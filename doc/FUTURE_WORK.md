# Future work

This document lists outstanding work that the paper-faithful bridge
implementation did NOT cover. Items are grouped by priority and
context.

## Priority 1 — must do before merging

### 1.1 Visual ParaView verification of bridge geometry

**What:** Open `out_postarchipelago.vtk` from
`./build/mesher_roi -p data/Agua.poly -u out -a 3 -T -J -K 2 -F 0.5 -L 1 -e -v`
in ParaView and verify visually that the 125 bridges:

1. Land physically between the two components they connect.
2. Do not extend outside the polyline.
3. Have correct orientation (not flipped).
4. Form a coherent chain when 2+ bridges connect the same components.

**Why this matters:** the quantitative metrics (549 cells, 110
components, centroid analysis) confirm the algorithm's correctness
but not the geometric placement. A bridge that connects via centroid
but extends 0.5×0.5 perpendicular to the components would still pass
the "is connected" check but produce a visually wrong mesh.

**Effort:** ~15 min on a workstation with ParaView installed.

**Owner:** the user (preferred — they have the actual visualization
tools and the eye for the paper's reference figures).

### 1.2 Re-run regression tests after final documentation commits

**What:** After all `doc/*.md` files are in their final state, run
`python3 scripts/test_bridge_clean_info.py build/mesher_roi` one more
time. Confirm 6/6 pass.

**Why this matters:** documentation changes only, but the final
documentation includes a `BUILD.md`-style section claiming "build is
clean". If a test starts failing, the claim is wrong.

**Effort:** ~30 s.

### 1.3 Build is clean confirmation

**What:** Run `cmake --build build` from a fresh checkout to confirm
no warnings or errors.

**Why this matters:** all 6 regression tests were last run when the
code was last rebuilt. Documentation changes don't affect the
binary, but a fresh build catches stale `.o` files.

**Effort:** ~30 s on the dev box.

## Priority 2 — important features

### 2.1 Command-line flag to disable bridge-joining (`-B`)

**What:** Add a CLI flag `-B` (default: false) that, when set to true,
skips the bridge-joining loop entirely. Useful for comparison:
running with and without `-J` on the same input should produce
visually different but quantitatively related results.

**Where:** `src/Mesher.cpp:118` (CLI parser) +
`src/Mesher.cpp:3392` (conditional skip around bridge loop).

**Effort:** ~30 min including the parser update.

**Why this matters:** the user wants to compare the old (pre-fix)
algorithm with the new one. Without `-B`, the only way to do this is
to revert the code, which is error-prone.

### 2.2 Override `preRefineForTusqh` base level

**What:** Add a CLI flag `-K0 <int>` (or similar) that overrides the
hardcoded `baseLevel=3` in `preRefineForTusqh`. Useful for users
whose initial grid is finer than 3.

**Where:** `src/Mesher.cpp:3622` (`baseLevel` parameter) +
`src/Mesher.cpp:118` (CLI parser).

**Effort:** ~20 min.

### 2.3 Pinch-template handling for 3D

**What:** Extend `bridgeSplitAtEdge` and the bridge loop to handle
3D pinch templates (paper §3.4 second half, briefly mentioned).

**Where:** currently 2D-only. The 3D path is in
`src/Quadrant.h:58-66` (where `Quadrant` has 8 corner points instead
of 4) but the bridge logic doesn't use them.

**Effort:** ~1 day.

**Why this matters:** the paper is 2D-only in this work, but the
codebase is set up for 3D. Future work needs to make the bridge
algorithm work in 3D too.

### 2.4 Speed up the BFS

**What:** The current BFS in `Mesher.cpp:3361-3520` is O(n) per
iteration but iterates the full `Quadrants` vector (which can have
~1400 entries for Chesapeake Bay at `-a 3`). For larger
`-a` values this could become slow.

A better approach: use `QuadTree`'s built-in BFS via the visitor
pattern instead of `unordered_map<unsigned int, unsigned int>`.

**Where:** `src/QuadTree.cpp:78-110` (existing
`visitEveryLeafQuadrant`) could be adapted.

**Effort:** ~4 hours.

### 2.5 Persistence diagrams (Fig 7 of paper)

**What:** Compute and emit persistence diagrams from the post-TUSQH
output, as in Fig 7 of the paper.

**Where:** new file `src/Persistence.cpp` (and add
`#include <gudhi/Persistence_representations.h>`).

**Dependencies:** Gudhi library (not currently in the project's
CMakeLists.txt). Would require adding it.

**Effort:** ~2-3 days.

**Why this matters:** Fig 7 is the paper's most striking
visualization and the user's main interest. But it requires an
external dependency (Gudhi) and substantial new code.

## Priority 3 — polish

### 3.1 Document `getExteriorDirection` and `isEdgeOnDomainBoundary` semantics

**What:** Add Doxygen comments to these two functions. Currently they
are described in `BUGS_FOUND.md` but not in the source.

**Where:** `src/Mesher.h:253-256` (declarations).

**Effort:** ~15 min.

### 3.2 Add a build flag `--enable-paper-faithful-bridge` (default on)

**What:** CMake option `MESHER_PAPER_FAITHFUL_BRIDGE` (default ON)
that gates the new code path. Default OFF would re-enable the
original isolated-quads algorithm (useful for A/B comparison).

**Where:** `CMakeLists.txt:21` + `src/Mesher.cpp:3361`.

**Effort:** ~20 min.

### 3.3 Replace `print` debugging with a proper logger

**What:** Currently `Mesher.cpp` has ~80 `cout` calls scattered
throughout. Replace with a `Logger` class that respects a `-v` flag
or similar.

**Where:** `src/Mesher.cpp` (every `cout << "    * ..."` line).

**Effort:** ~3 hours.

### 3.4 Document the reclassification step

**What:** The bridge quads come out as `WindingState::Unknown` after
`bridgeSplitAtEdge`. They are reclassified via
`WindingNumberVisitor` afterward (`Mesher.cpp:3560-3590`). This is
non-obvious; document in TUSQH_WINDING_GUIDE §3.11.9.

**Effort:** ~10 min.

## Priority 4 — known issues

### 4.1 `Quadrant::getIndex()` const-ness (Issue #7 in BUGS_FOUND)

**What:** Declare as `const unsigned int &getIndex() const;` so it
works on `const Quadrant&`.

**Where:** `src/Quadrant.h:133`.

**Effort:** ~5 min.

**Caveat:** changing it might break callers that don't expect a
const-qualified method. Audit all callers first.

### 4.2 `generateQuadtreeMesh` argument order

**What:** Currently `Mesher.cpp:271` calls with arguments in
wrong order. Only triggers in legacy non-TUSQH path.

**Where:** `src/Mesher.cpp:271`.

**Effort:** ~10 min.

**Caveat:** fixing this would change legacy output meshes. May be
intentional; investigate before changing.

### 4.3 Committed LaTeX build artifacts

**What:** `features/roi2D.bbl` and `features/roi2D.blg` should be in
`.gitignore`, not committed.

**Where:** `features/.gitignore` (does not exist).

**Effort:** ~5 min.

### 4.4 CMake warnings about deprecated `find_package` calls

**What:** Several `find_package(X REQUIRED)` calls in `CMakeLists.txt`
emit "X not found, trying without REQUIRED" warnings on clean
checkouts. Should be cleaned up.

**Where:** `CMakeLists.txt:32-58`.

**Effort:** ~30 min.

## Priority 5 — research

### 5.1 Adaptive thresholding

**What:** The current threshold `-F` is uniform across the entire
mesh. The paper hints at adaptive thresholds in §3.4 (briefly
mentioned). Implementing this would produce meshes that better
preserve local topology.

**Where:** `Mesher::windingSubdivide` and `preRefineForTusqh`.

**Effort:** ~1 week.

### 5.2 Multi-resolution output

**What:** Currently the output is one mesh. The paper discusses
hierarchical output (different resolutions for different regions).

**Where:** new code, possibly `src/Hierarchy.cpp`.

**Effort:** ~2 weeks.

### 5.3 Comparison with Gudhi-based TUSQH

**What:** Compare our TUSQH implementation with the original
Gudhi-based reference (paper §4, reference impl on GitHub).

**Where:** new benchmark `benchmarks/tusqh_compare.cpp`.

**Effort:** ~1 week (mostly setting up the benchmark infrastructure).

## Done — already implemented

For completeness, here is what the branch ALREADY has:

1. ✅ Sub-cell volume fractions (paper §3.3).
2. ✅ Full-grid BFS for bridge-joining (paper §3.4 first part).
3. ✅ Bridge templates that connect components (paper §3.4 second part).
4. ✅ Domain boundary exclusion (Issue #1).
5. ✅ Correct outward direction (Issue #2).
6. ✅ MapEdges cleanup (Issue #3).
7. ✅ AllOutside preservation (Step 1).
8. ✅ Output filtering (Step 4).
9. ✅ Pre-refinement for coarse grids (Issue #6).
10. ✅ Empty-Quadrants guard (Issue #5).
11. ✅ Quantitative validation script.
12. ✅ 6 regression tests.
13. ✅ Documentation overhaul (this and 7 STEP_*.md files).
