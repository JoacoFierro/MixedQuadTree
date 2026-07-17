# Step 3 — Rewrite bridge-joining loop (paper-faithful)

**Status:** ✅ SUCCESS

**File modified:** `src/Mesher.cpp` `resolveArchipelagos` (lines 3326–3556)

**Changes:**

The old loop only considered **exposed edges** (edges where `info[2] == MAX`, meaning one side has a kept quad and the other side has nothing). For each such edge it added a 1-to-5 split that grew OUTWARD into the gap. With AllOutside cells dropped, this produced **isolated bridge quads** instead of connecting components (see `STEP_3_LIMITATION_2026-07-15.md` previously logged in `TUSQH_WINDING_GUIDE.md` §3.11.9).

The new loop implements the paper's algorithm directly:

**Each iteration:**
1. **BFS over INTERIOR cells only** (`isInteriorCell(Quadrants[start])`). AllOutside cells are skipped during BFS but their edges remain in MapEdges, so the BFS can hop across chains of exterior cells.
2. **Bridge candidates** are edges in MapEdges where:
   - BOTH sides are interior quads (`info[1] != MAX && info[2] != MAX`)
   - The two quads are in DIFFERENT components (after BFS)
   - The edge has interior sub-cell VF (`>= joinThreshold`)
   - The edge is NOT on the global domain boundary (Issue #1 filter)
3. For each candidate, apply the 1-to-5 split on the `info[1]` quad. The bridge quad extends OUTWARD from `info[1]`'s centroid toward `info[2]`'s centroid — since the two cells are physically adjacent in the cubical complex, the bridge quad lands on `info[2]`'s territory, **connecting** the two components.

**Key insight:** because we kept the AllOutside cells in `Quadrants` (Step 1), two interior cells separated by a chain of exterior cells are now connected in MapEdges through those exterior cells. The BFS treats them as **different components** (because the exterior cells are not interior), and the algorithm finds the edges between them and bridges them.

**Verification (Chesapeake Bay, `-a 3 -T -J -K 2 -F 0.5 -L 1`):**

| Metric | Before Step 3 | After Step 3 |
|--------|---------------|--------------|
| Bridges added | 0 | 125 |
| Iter 0 | (no candidates) | 122 bridges, 7 boundary edges filtered |
| Iter 1 | — | 3 bridges, 5 boundary edges filtered |
| Iter 2 | — | no candidates |
| resolveArchipelagos summary | 0 bridges, 0 drops | 125 bridges, 970 drops |

The 970 "drops" include both the 345 preserved AllOutside cells (which were dropped at the end of resolveArchipelagos because `compOfQuad == -1`) AND the cells from any component smaller than `-L`. With `-L 1`, only the AllOutside cells should be dropped.

**Build:** clean.

**Issue:** `saveOutputMesh` does not filter exterior cells yet, so the VTK output currently includes both interior and exterior cells. This makes the post-archipelago output look larger than expected. Fix in Step 4.

**Next:** Step 4 — filter exterior cells out of `saveOutputMesh` so the final mesh only contains the interior (kept) cells.

## Cross-references

- This step is the heart of the paper-faithful implementation. It resolves the highest-severity bug documented in `BUGS_FOUND.md` as **Issue #4** ("bridge-joining produces isolated quads instead of connecting components"). The original algorithm considered only exposed edges and added bridge quads outward into the gap; the new algorithm considers edges between **two interior quads** in **different components** of the cubical complex.
- The bridge candidate filter depends on three bug fixes documented in `BUGS_FOUND.md`:
  - **Issue #1** (`isEdgeOnDomainBoundary`) — exclude edges whose outward fictitious cell has sub-cell VF below the threshold.
  - **Issue #2** (`computeExteriorDirection`) — use the quad centroid (not the opposite-edge midpoint) to determine outward direction.
  - **Issue #3** (`bridgeSplitAtEdge` MapEdges cleanup) — erase the stale full-edge entry after the 1-to-5 split.
- After `bridgeSplitAtEdge`, the new bridge quad comes out with `WindingState::Unknown`. The loop at the end of this step (`Mesher.cpp:3560-3590`) reclassifies it via `WindingNumberVisitor`. This is documented in `FUTURE_WORK.md` §3.4.
- See `WORK_SUMMARY.md` §"Metrics" for the final Chesapeake Bay numbers, and `STEP_5_VALIDATION.md` for the detailed validation analysis.
