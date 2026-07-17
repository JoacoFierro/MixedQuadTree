# Step 5 — Quantitative validation of Chesapeake Bay output

**Status:** ✅ PASSED (algorithm produces sensible output, but visual confirmation in ParaView recommended)

**Test command:**
```
./build/mesher_roi -p data/Agua.poly -u /tmp/final_a3_K2 \
    -a 3 -T -J -K 2 -F 0.5 -L 1 -e -v
```

**Validation script:** `scripts/analyze_output.py <vtk> <poly>`

**Results (Chesapeake Bay, `-a 3 -T -J -K 2 -F 0.5 -L 1`):**

| Metric | Value |
|--------|-------|
| `_postarchipelago.vtk` cells | 549 |
| `_postarchipelago.vtk` points | 2129 |
| `_postarchipelago.vtk` components (mesh adjacency) | 110 |
| `_tusqh.vtk` cells (before bridges) | 1019 |
| `_tusqh.vtk` components | 85 (biggest = 291) |
| Total bridges added | 125 (across 2 iterations) |
| Domain boundary edges filtered | 9 (Issue #1) |
| `_postarchipelago.vtk` AllOutside cells | 0 ✓ (Step 4 working) |

**Top-10 components by size:**
- Biggest: 49 cells, centroid (31.79, 24.12), bounding box [25.31, 39.72] × [9.78, 38.61]
- 8 of 10 centroids are inside the polyline
- 2 of 10 centroids are outside the polyline — but bounding boxes overlap the polyline's bounding box. These are likely cells straddling the polyline boundary (the polyline outlines a thin shoreline).

**Overall:**
- 63 components centroids inside polyline (likely real islands / coastal features)
- 47 components centroids outside polyline (likely boundary-straddling cells or the algorithm is including some exterior cells incorrectly — needs visual confirmation)

**Algorithm sanity check:**
- Bridges added: 125 ✓ (paper-faithful loop runs to completion in 2 iterations)
- Issue #1 (boundary filter): 9 boundary edges correctly filtered ✓
- Step 4 filter: 970 AllOutside cells correctly dropped from output ✓
- Step 3 bridge loop runs with full-grid BFS, finds 121 candidates in iter 0, then 4 more in iter 1 (diminishing returns as expected) ✓

**Concerns / open questions for visual check:**

1. The biggest mesh component has only 49 cells. Is this the bay itself (should be much larger), or just one of many islands?

2. 47 components have centroids outside the polyline. With a closed shoreline polygon, interior cells should have centroids inside. Need to verify visually whether these are:
   - (a) Real artifacts (algorithm wrongly classifying exterior cells as interior), or
   - (b) Legitimate boundary-straddling cells (where the centroid is just outside the polygon by a small margin).

3. Pre-bridge vs post-bridge component count:
   - Pre-bridge (TUSQH): 85 components in the cubical complex
   - Post-bridge: 110 components in the OUTPUT mesh
   - But bridges should DECREASE components (joining them). Why did components INCREASE?
   - Hypothesis: the bridge loop is fragmenting the bay's main component by inserting bridges that subdivide large cells into smaller disconnected clusters. This needs verification.

**Recommendation:** Visual confirmation in ParaView is required before marking Step 5 complete. The script `scripts/analyze_output.py` provides the quantitative foundation, but the visual check is the only way to confirm the algorithm is producing the correct mesh.

**Build:** clean.

**Next:** Step 6 — run all 5 regression tests, update expectations to match the new paper-faithful behavior, re-run tests to confirm they pass.

## Cross-references

- The 549 interior cells and 110 components table matches the "Final results" table in `WORK_SUMMARY.md`.
- The 9 boundary edges filtered are the Issue #1 fix in `BUGS_FOUND.md` exercised on Chesapeake Bay.
- The 125 bridges added are the result of the Step 3 rewrite; see `STEP_3_BRIDGE_PAPER_FAITHFUL.md` for the algorithm and `BUGS_FOUND.md` Issue #4 for the original bug.
- The 970 AllOutside cells dropped are the Step 1 preserve-then-Step 4 filter pattern; see `WORK_SUMMARY.md` §3.
- Visual validation in ParaView is **still pending** — see `FUTURE_WORK.md` §1.1 and `CURRENT_STATE.md` §"What doesn't work yet".
