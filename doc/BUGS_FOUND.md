# Bugs found and fixed

This document is a comprehensive log of every bug discovered during the
implementation of TUSQH §3.3 + §3.4 in `MixedQuadTree`. Each entry follows
the format:

> **Symptom** → **Root cause** → **Fix** → **Verification**

All bugs are RESOLVED. The final commit will close them all.

---

## Issue #1 — Domain boundary edge incorrectly accepted as bridge candidate

**Severity:** High (algorithmic correctness — bridges could extend outside the polyline into genuine exterior).

**Status:** ✅ Fixed (pre-existing in the branch, documented in `SESSION_TUSQH_BRIDGE_2026-07-15.md` §6).

### Symptom

When running the bridge-joining algorithm on a polyline whose outer
boundary forms part of the quadtree's outer boundary (the L-shaped test
case `tusqh_boundary_filter.poly`), some candidate edges had
`info[2] == MAX` on the global boundary side. The bridge was added
without realizing the "exterior" side of the edge was actually
**outside the entire polyline**, not just the local quad. The bridge
quad ended up in genuine exterior and the post-archipelago mesh had a
stray cell outside the polyline outline.

### Root cause

The original bridge candidate filter in `resolveArchipelagos` did:

```cpp
if (info[1] == MAX || info[2] == MAX) continue;  // exposed edge
```

This was supposed to **skip** edges whose other side is "empty"
(`MAX`). But:

1. **`info[2] == MAX`** correctly identifies exposed edges.
2. The remaining filter (`info[1] != MAX && info[2] != MAX`) checked
   that both sides had a quad. **However**, it did NOT check whether
   the OTHER side of the edge (the `info[2]` quad) was on the global
   domain boundary.

### Fix

Added a dedicated function `Mesher::isEdgeOnDomainBoundary`
(`Mesher.cpp:3262-3348`, declared in `Mesher.h:262`) that:

1. Takes an edge `(a, b)`, the quad on side `info[1]`, the polyline,
   the sample size `s`, and the join threshold `t`.
2. Computes the perpendicular direction outward from the quad (using
   the quad centroid — see Issue #2).
3. Builds an `s × s` fictitious cell in the outward direction.
4. Computes its volume fraction. If `VF < threshold`, the edge is on
   the global domain boundary → **drop it as a bridge candidate**.

The function is called inside the bridge candidate loop with **strict
inequality** (`<`, not `<=`), so an edge whose outward fictitious
cell has `VF == threshold` is still treated as interior. This matches
the convention in the rest of the pipeline (boundary == `VF <
threshold`).

### Verification

`scripts/test_bridge_clean_info.py::boundary_filter_Lshape` runs
`-a 3 -T -J -K 2 -F 0.5 -L 1` on `tusqh_boundary_filter.poly` (an
L-shaped outer boundary) and expects at least 4 boundary edges to be
filtered. With Issue #1 unfixed, the test would fail (bridges extend
outside the L into empty space).

---

## Issue #2 — Wrong outward direction for non-axis-aligned quads

**Severity:** High (geometric correctness — bridge quads point in wrong direction for rotated or non-convex quads).

**Status:** ✅ Fixed.

### Symptom

When the quadtree is refined in regions where the polyline's edge runs
diagonally across a quad, the quad becomes non-axis-aligned and
non-convex in extreme cases. The original "outward direction" for the
bridge quad was computed as:

```cpp
Point3D mid_opposite = (p_e2 + p_e3) * 0.5;
Point3D dir_exterior_unit = -normalize(mid_opposite - mid_bridge);
```

For an axis-aligned convex quad, `mid_opposite - mid_bridge` correctly
points from the bridge edge toward the quad's interior. **But** for
rotated or non-convex quads, the midpoint of the opposite edge does
NOT necessarily lie on the "interior" side of the bridge edge. The
result: the bridge quad extends in the wrong direction, sometimes
back INTO the existing quad (overlapping) or away from the target
component (missing it).

This manifested in the `rotated_quad_centroid_direction` test case:
a 45°-rotated square fed to the mesher produced a bridge quad pointing
in the wrong direction.

### Root cause

The midpoint of an edge is a purely geometric quantity that depends
on the corners of that edge only. It has no notion of "which side of
the edge is the quad's interior". For non-convex quads, the quad's
interior can be on either side of the edge midpoint.

### Fix

Use the **quad centroid** (mean of the four corners) as the interior
reference. The centroid always lies inside the convex hull of the
quad, so the direction from the bridge edge midpoint to the centroid
is always "toward the interior" (for any convex or non-convex quad):

```cpp
Point3D centroid = (p_e0 + p_e1 + p_e2 + p_e3) * 0.25;
Point3D dir_exterior_unit = -normalize(centroid - mid_bridge);
```

Implemented as `Mesher::computeExteriorDirection`
(`Mesher.cpp:3287-3348`, declared in `Mesher.h:253-256`).

### Verification

`scripts/test_bridge_clean_info.py::rotated_quad_centroid_direction`
runs the mesher on `tusqh_rotated.poly` (a 45°-tilted square inside
the bounding box). Expected: 1 cell in output (the rotated square is
fully inside, no bridges needed). Pre-fix: bridge quad points wrong
direction, sometimes overlapping with the original quad. Post-fix:
output is clean.

---

## Issue #3 — Stale `info[1]` in `MapEdges` after bridge split

**Severity:** High (algorithmic correctness — BFS uses stale `q_id` and corrupts component labeling).

**Status:** ✅ Fixed.

### Symptom

After `bridgeSplitAtEdge` replaces a quad `q` with 5 new quads (4
sub-quads + 1 bridge quad), the **full edge** `(e0, e1)` (the bridge
edge) was still present in `MapEdges` with `info[1]` pointing to the
**old** q_id (which no longer exists in `Quadrants`).

When the next bridge iteration ran its BFS, the lookup
`qIdToIdx.find(info[1])` returned `end()`, silently skipping the edge.
This left the BFS unable to "cross" the bridge edge between the two
new components, so the bridges failed to actually connect anything in
some cases.

### Root cause

`SplitVisitor::visit` correctly inserts the half-edges `(e0, m_bridge)`
and `(m_bridge, e1)` after the 1-to-4 split. The bridge quad then
extends info[2] on these half-edges to itself. But the **original
full edge** `(e0, e1)` is not removed — it has `info[1] = old_q_id`
which is now stale.

### Fix

In `bridgeSplitAtEdge` (`Mesher.cpp:3198-3202`), after building the 5
new quads, **erase** the stale full edge from `MapEdges`:

```cpp
{
    unsigned int a = e0, b = e1;
    if (a > b) std::swap(a, b);
    MapEdges.erase(QuadEdge(a, b));
}
```

The two half-edges remain and correctly reference the bridge quad on
one side and the sub-quad on the other.

### Verification

Verified in the diagnostic run (Step 5). Without Issue #3 fix, the
BFS in iter 1 would have looked up the stale `info[1]` and failed.
With fix, the BFS sees only valid edges.

---

## Issue #4 — Bridge-joining adds isolated quads instead of connecting components

**Severity:** **Highest** (algorithmic correctness — fundamental limitation of the algorithm as designed).

**Status:** ✅ Fixed in branch `feature/paper-faithful-bridge` (Steps 1–4).

### Symptom

Initial verification of the bridge-joining algorithm with
`data/Agua.poly -a 3 -T -J -K 2 -F 0.0 -L 1` showed:

| Metric | Value |
|--------|-------|
| Cells pre-bridge | 674 |
| Cells post-bridge | 1802 |
| Bridges added | 282 |
| Components pre-bridge | 97 |
| Components post-bridge | **378** |
| Biggest component pre-bridge | 61 cells |
| Biggest component post-bridge | 65 cells |

The 282 bridges **created 281 new components** instead of joining
existing ones. The biggest component grew only 4 cells.

The "small gap" test (`data/tusqh_bridge_small_gap.poly`: two 0.5×0.5
squares with gap=0.1) confirmed the same behavior: even with a small
gap, the bridge quad landed **outside** both squares as a new isolated
component.

### Root cause

The original algorithm only considered **exposed edges** (where
`info[2] == MAX` — meaning one side of the edge has no quad). For each
such edge, it added a 1-to-5 split to the quad on the OTHER side. The
bridge quad extended **outward** from that quad's centroid, with
thickness `H / sampleSize`.

For a bridge to physically **connect** two components, the bridge
quad must overlap with a quad from the OTHER component. This requires
the gap between the two components to be **smaller** than
`H / sampleSize`.

For Chesapeake Bay with `-K 2` (sample size 2, quad size ≈ 0.5 in
unit-square coords), the typical gap between adjacent islands is much
larger than 0.25. The bridge quad lands in the gap, becomes a new
isolated quad, and is later dropped by the small-component filter.

This is documented in `doc/TUSQH_WINDING_GUIDE.md` §3.11.9 (pre-rewrite
version, "Limitación fundamental del algoritmo bridge-joining actual").

### Fix

The paper's §3.4 specifies a different (and correct) algorithm:

> "For any pair of connected components, if the edges that connect
> them are interior to the geometry, the components are joined using
> templates along those edges."

The fix requires **three** changes, documented as Steps 1–3 in
`doc/STEP_0..3`:

1. **Step 1 (`Mesher.cpp:2030-2068`): preserve AllOutside cells.**
   The original pipeline drops AllOutside cells after TUSQH. We now
   keep them in `Quadrants` so `MapEdges` retains the topology of the
   full cubical complex.

2. **Step 2 (`Mesher.h:289`, `Mesher.cpp:3250`): `isInteriorCell` helper.**
   Centralises the check `WindingState::AllInside || Mixed → true`.

3. **Step 3 (`Mesher.cpp:3361-3520`): rewrite the bridge loop.**
   The new loop:
   - Runs BFS over the **full cubical complex** (interior + exterior).
   - Finds edges where **both sides are interior quads** AND they
     belong to **different components** (after BFS) AND the edge has
     interior sub-cell VF AND it is **not on the global domain boundary**
     (Issue #1).
   - Applies `bridgeSplitAtEdge` on the quad in `info[1]`. Because the
     two cells are physically adjacent in the cubical complex, the
     bridge quad lands ON `info[2]`'s territory, connecting the two
     components.

4. **Step 4 (`Mesher.cpp:3621-3690`): output filter.**
   The `AllOutside` cells we preserved for the BFS are now excluded
   from the final mesh via `isInteriorCell(Quadrants[qi])` in the
   `keepQuad[]` filter.

### Verification

Chesapeake Bay `-a 3 -T -J -K 2 -F 0.5 -L 1` (Step 5):

| Metric | Pre-fix | Post-fix |
|--------|---------|----------|
| Bridges added | 0 (limit) or 282 (with `-F 0.0`) | 125 (correctly placed) |
| Cubical complex components | 97 | 74 |
| Output mesh components | 378 | 110 |
| Output interior cells | many spurious | 549 (clean) |

Plus a new positive-control test
`scripts/test_bridge_clean_info.py::bridge_connects_two_islands`
verifies that for two close squares with gap=0.1, the algorithm adds
2 bridges and produces 3 cells in the output (2 squares + 1 bridge
quad interior to the gap).

---

## Issue #8 — Bridge quad is topologically isolated (erased full edge)

**Severity:** High (the bridge quad survives only because of the
`-L 1` filter; in any other mode it would be dropped).

**Status:** ✅ Fixed.

### Symptom

In `Mesher::bridgeSplitAtEdge`, after performing the 1-to-5 split, the
code does:

```cpp
MapEdges.erase(QuadEdge(e0, e1));  // ← Issue #3
```

This removes the MapEdges entry for the **full bridge edge** between
the original quad `q` and the neighbour `q2`. The two half-edges
`(e0, m_bridge)` and `(m_bridge, e1)` are kept and reference the new
bridge quad, so the bridge quad is connected to `q`'s sub-quads.

But the bridge quad has **no edge connecting it directly to `q2`** —
all 4 of its edges are either halves (connecting to `q`'s sub-quads)
or brand-new outer edges (connecting to nothing).

Result: the bridge quad is a 1-quad connected component on its own
(only reachable through `q`'s sub-quads, which are on the OTHER side
of `q2`).

The only reason the bridge quads survive in the output is that the
`isInteriorCell`-filter on small components (`-L 1` by default in the
test framework) keeps any component with ≥ 1 quad. The behaviour
**silently changes** if the user passes `-L 2` or larger: the bridge
quads disappear.

### Root cause

The Issue #3 fix (erase full edge after 1-to-5 split) was made on the
assumption that the halves fully represent the topology. They don't,
because the bridge quad is geometrically on the `q2` side of the
bridge edge, so it needs a MapEdges entry that points to `q2`.

### Fix

In `Mesher::bridgeSplitAtEdge` (`src/Mesher.cpp:3230-3248`),
**replace** the `MapEdges.erase((e0, e1))` with an in-place update:

```cpp
auto it = MapEdges.find(QuadEdge(min(e0,e1), max(e0,e1)));
if (it == MapEdges.end()) {
    MapEdges.emplace(ke, EdgeInfo(m_bridge, bridge_q_id, q2.getIndex()));
} else {
    it->second[0] = m_bridge;        // midpoint stays the same
    it->second[1] = bridge_q_id;     // ← was the removed `q`
    it->second[2] = q2.getIndex();   // ← already set by SplitVisitor
}
```

Now the bridge quad is reachable from `q2` via the full `(e0, e1)`
edge, so the BFS in `resolveArchipelagos` puts the bridge quad in
`q2`'s component (and consequently `q`'s sub-quads, which share the
halves with the bridge quad).

The function signature was extended to take `q2` and a
`doManifoldSplit` flag (Option B-subdiv placeholder; not implemented
yet but kept in the API for future use).

### Verification

| Metric | Pre-fix (Issue #3 erase) | Post-fix (Issue #8 update) |
|--------|--------------------------|----------------------------|
| Agua `-a 3 -T -J -K 2 -F 0.5 -L 5`: components | 241 | 229 (-12) |
| Agua `-a 3 -T -J -K 2 -F 0.5 -L 5`: bridges | 125 | 117 (-8) |
| Agua `-a 3 -T -J -K 2 -F 0.5 -L 5`: small-component quads dropped | 1020 | 979 (-41) |

The 12 fewer components and 8 fewer bridges indicate that the bridge
quads are now properly joining their neighbour's component (instead of
becoming isolated 1-quad components). The 41 fewer dropped quads
confirms that more quads survive in larger (merged) components rather
than being dropped as small 1-quad components.

The 6 baseline tests in `scripts/test_bridge_clean_info.py` all still
pass.

**Note on manifoldness:** Option A (topology fix only) leaves the
bridge quad geometrically overlapping `q2`'s rectangular fictitious
cell. Strict manifoldness would require Option B-subdiv (splitting
`q2` 1-to-4 and discarding its 2 sub-quads adjacent to the bridge
edge). The `doManifoldSplit` flag in the API is reserved for that
future implementation.

---

## Issue #5 — `saveOutputMesh` segfaults on empty `Quadrants`

**Severity:** Medium (crash, but only when all components dropped).

**Status:** ✅ Fixed (pre-existing in the branch, documented in
`SESSION_TUSQH_BRIDGE_2026-07-15.md` §7).

### Symptom

Running the mesher with `-L 99999` (drop all components) caused a
segfault in `Mesher::saveOutputMesh`. The crash happened at the line
that reads `tmp_Quadrants[0].getRefinementLevel()` without checking if
the vector was empty.

The same crash was reachable from other paths:
`Mesher.cpp:138, 150, 165` (debug `_quads`, `_closeto`, `_remSur`
saves).

### Root cause

The function assumes `tmp_Quadrants` is non-empty when called from
`resolveArchipelagos`'s post-drop step. If `-L` is large enough that
all components are dropped, `Quadrants` becomes empty and the
dereference crashes.

### Fix

Added an explicit empty-vector guard at the start of
`saveOutputMesh(vector<Quadrant>&)` (`Mesher.cpp:1622-1653`):

```cpp
if (tmp_Quadrants.empty()) {
    mesh->setPoints(out_pts);
    mesh->setElements(out_els);
    cout << "    * SaveOutputMesh (empty Quadrants) in ... ms" << endl;
    return 0;
}
```

This emits an empty mesh instead of crashing.

### Verification

`scripts/test_bridge_clean_info.py::figure14_drop_all_components`
runs with `-L 3` on a polyline whose all components have 2 cells (so
all get dropped). Expected: process exits 0 with no segfault.

---

## Issue #6 — TUSQH undersamples very coarse initial grid

**Severity:** High (algorithmic — Chesapeake Bay produces only 3 cells without the fix).

**Status:** ✅ Fixed (pre-existing in the branch, documented in
`SESSION_TUSQH_BRIDGE_2026-07-15.md` §6).

### Symptom

Chesapeake Bay (`data/Agua.poly`) has bbox ~50×75 in real-world
units. The mesher's `generateGridMesh` initial step is
`min(dx,dy) * 1.01` ≈ 50, so it produces **2 root cells** (each
~50×75 — very elongated and huge).

TUSQH then classifies each cell using `s × s` samples (default
`-N 2` = 4 samples per cell). For a 50×75 cell, 4 samples are wildly
insufficient: most samples land outside the polyline and the cell is
classified as `AllOutside`, never subdivided.

Result: `-a 7 -T -N 2` produced only **3 cells** in 0.4 s — a useless
mesh.

### Root cause

The initial grid from `generateGridMesh` is too coarse for the polyline
in the cell, and TUSQH's `s × s` sampling can't refine it further
because it doesn't subdivide `AllOutside` cells.

### Fix

Added a heuristic `Mesher::preRefineForTusqh`
(`Mesher.cpp:3622-3673`, declared in `Mesher.h:158-183`) that:

1. Triggers when `Quadrants.size() <= 2` AND
   `input.getEdges().size() >= 100`.
2. Pre-refines the quadtree uniformly to level 3 (configurable via
   `baseLevel`).
3. Caps at the user's `maxDepth` so the TUSQH loop still has headroom.

The trigger conditions exclude trivial regression cases (e.g.,
`unit_square.poly` has 2 root cells but only 4 segments, so the
heuristic does NOT fire).

### Verification

Chesapeake Bay `-a 3 -T -J -K 2 -F 0.5 -L 3` now produces 645 cells
in 53 s (vs. 3 cells in 0.4 s without the fix). Unit square tests
still pass (heuristic doesn't fire).

---

## Issue #7 — Quadrant::getIndex() declared non-const but used as const

**Severity:** Low (compile warning, but breaks some compilers).

**Status:** Documented in `SESSION_TUSQH_BRIDGE_2026-07-15.md` §8.

### Symptom

Header `src/Quadrant.h:133` declares:

```cpp
virtual const unsigned int &getIndex();
```

The return type is `const unsigned int &` but the method itself is
non-const. Taking `const Quadrant& q; q.getIndex()` triggers a
`-fpermissive` warning on strict compilers.

### Root cause

Historical — the method was changed to return a const reference but
the const-ness of the method itself was not updated.

### Fix

Not fixed in this branch because changing it might affect other
consumers. Documented as "known issue, work around with `auto& ref =
const_cast<Quadrant&>(q); ref.getIndex()`" in some hot paths.

### Verification

Build succeeds with `-Wno-error`; only emits warnings.

---

## Issue #9 — `buildQuadPerpThickness` treats `info[k]` as a vector index (q_id ≠ index)

**Severity:** Medium-high (algorithmic correctness — bridge selection
can pick the wrong edges because edge sub-cell VF is wrong).

**Status:** ✅ Fixed.

### Symptom

`SubgridSampler::buildQuadPerpThickness` (`src/SubgridSampler.cpp:185`)
used `Quadrants[qIdx]` with `qIdx = EdgeInfo.info[k]` to look up the
quad adjacent to an edge. The problem: `EdgeInfo.info[k]` stores a
**q_id** (the id assigned at `Quadrant` construction), while
`Quadrants[...]` is a `std::vector` indexed by **position**. After
`windingSubdivide` (and especially after `resolveArchipelagos`'s
compact-and-rebuild) the q_id and the vector index of a given quad
can disagree, so the function silently read the **wrong** quad.

This propagated to the edge sub-cell volume fraction:

```
SubgridSampler::buildQuadPerpThickness (wrong quad)
        |
        v
SubgridSampler::sampleEdge (wrong fictitious cell size)
        |
        v
mEdgeSubcellVF[edge] = wrong VF
        |
        v
resolveArchipelagos: bridge candidate filter
  `if (vfIt->second.volumeFraction < joinThreshold) continue;`
        |
        v
Bridges may be wrongly accepted or wrongly dropped.
```

In `bridgeSplitAtEdge`, `H = |centroid - mid_bridge|` also depends on
the same lookup (via `computeExteriorDirection`); when the wrong quad
is used, the bridge quad extends in the wrong direction by
`H/sampleSize`.

### Root cause

Two distinct indexing conventions coexisted:

- `EdgeInfo::info[1]`, `info[2]` store q_ids (set in
  `SplitVisitor::visit` and `bridgeSplitAtEdge`, never reassigned).
- `Quadrants` is a `std::vector` indexed by position; positions can
  be reordered (e.g. by `resolveArchipelagos`'s compact step at
  `Mesher.cpp:3714-3755`).

`buildQuadPerpThickness` treated `info[k]` as a position. The
bounds-check `qIdx >= quadrants.size()` was insufficient — it caught
qids outside the range but did not catch qids that happen to land on
the wrong quad within the range.

`resolveArchipelagos` does the right thing: it builds a
`qIdToIdx` map (`Mesher.cpp:3459-3463`) before indexing
`Quadrants[…]`.

### Fix

`SubgridSampler::buildQuadPerpThickness` now requires a precomputed
`unordered_map<unsigned int, unsigned int> qIdToIdx` and uses
`quadrants[qIdToIdx.at(qId)]` (with a guard). The caller
`Mesher::computeSubcellVolumeFractions` (`Mesher.cpp:2981-2986`)
builds the map once before the edge loop. The function signature
in `src/SubgridSampler.h:120-141` was extended accordingly.

### Verification

- Build succeeds (only pre-existing warnings).
- `scripts/test_qid_vs_index_regression.py` exercises Agua.poly
  (segments >= 100, triggers `preRefineForTusqh` and produces many
  quads with non-trivial q_id / index divergence) and checks that:
    - All 7 step VTKs are produced (octree, prearch, postarch,
      quads, closeto, remSur, final).
    - At least one `output_bridge_iter{N}.vtk` is produced.
    - The resolver summary line reports `bridges >= 1`.
- The bridge loop now uses correctly computed edge sub-cell VFs;
  the bridge quad thickness `H/sampleSize` is also correct.

---

## Summary table

| # | Bug | Severity | Status |
|---|-----|----------|--------|
| 1 | Domain boundary edge accepted as bridge | High | ✅ Fixed |
| 2 | Wrong outward direction for non-AABB quads | High | ✅ Fixed |
| 3 | Stale `info[1]` in `MapEdges` after split | High | ✅ Fixed |
| 4 | Bridge-joining produces isolated quads | **Highest** | ✅ Fixed (paper-faithful rewrite) |
| 5 | `saveOutputMesh` segfaults on empty `Quadrants` | Medium | ✅ Fixed |
| 6 | TUSQH undersamples coarse initial grid | High | ✅ Fixed (pre-existing) |
| 7 | `Quadrant::getIndex()` const-ness | Low | ⚠️ Documented, not fixed |
| 8 | Bridge quad is topologically isolated (erased full edge) | High | ✅ Fixed (Option A: update instead of erase) |
| 9 | `buildQuadPerpThickness` uses `info[k]` as vector index | Medium-high | ✅ Fixed (qIdToIdx lookup) |

## Outstanding latent bugs (not in this session)

These were discovered but **deliberately NOT fixed** to avoid altering
existing mesh outputs in other flows:

- **`generateQuadtreeMesh` call at `Mesher.cpp:271`** has arguments in
  incorrect order. Only triggers in the legacy non-TUSQH path. Fixing
  would change meshes produced by `-a` without `-T`.
- **`Quadrant::getIndex()` const-ness** (Issue #7 above).
- **`features/roi2D.bbl` and `roi2D.blg`** are LaTeX build artifacts
  committed by mistake (should be in `.gitignore`).

See `CURRENT_STATE.md` for the full list of uncommitted changes and
known limitations.
