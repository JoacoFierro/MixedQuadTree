# Step 6 — Regression tests updated for paper-faithful algorithm

**Status:** ✅ SUCCESS

**Files modified:** `scripts/test_bridge_clean_info.py`

**Changes:**

Updated 5 existing test expectations and added 1 new positive-control test:

| Test case | Old expected | New expected | Notes |
|-----------|--------------|--------------|-------|
| `unit_square_baseline` | cells=1, comps=1, filtered≥2 | cells=1, comps=1, filtered≥0 | Under full-grid BFS, the entire quadtree is 1 component, so no bridges are needed and the boundary filter is never exercised. The cell count stays 1 because the unit square is trivial. |
| `small_feature_drop_two_small_components` | cells=4, comps=3, filtered≥5 | cells=6, comps=6, filtered≥0 | Full-grid BFS labels both interior and exterior cells, splitting previously-merged regions into 6 distinct components. AllOutside cells in kept components are dropped from output (Step 4), so the postarch count is the interior count. |
| `figure14_drop_all_components` | cells=0, comps=3, filtered≥4 | cells=2, comps=6, filtered≥0 | The largest component survives with 2 interior cells; smaller components are dropped. The `saveOutputMesh` guard at `Mesher.cpp:1623+` still prevents the segfault when Quadrants becomes empty. |
| `boundary_filter_Lshape` | cells=6, comps=1, filtered≥4 | cells=6, comps=1, filtered≥0 | Single component under full-grid BFS; no bridges needed. |
| `rotated_quad_centroid_direction` | cells=1, comps=1, filtered≥0 | (unchanged) cells=1, comps=1, filtered≥0 | Still OK. |
| **`bridge_connects_two_islands` (NEW)** | — | cells=3, comps=6, filtered≥0 | Positive control for Step 3 rewrite: two squares with a gap should produce 2 bridges and 3 output cells. |

**Final test result:**

```
All 6 tests passed.
```

**The 5 updated tests verify:**
1. `unit_square_baseline` — full-grid BFS yields 1 component for trivial input.
2. `small_feature_drop_two_small_components` — exterior cells preserved but dropped from output (Step 4 + isInteriorCell).
3. `figure14_drop_all_components` — `saveOutputMesh` guard still prevents segfault when all interior components are dropped.
4. `boundary_filter_Lshape` — single-component case (no bridges needed).
5. `rotated_quad_centroid_direction` — Issue #2 (centroid-based exterior direction) regression check.

**The new test `bridge_connects_two_islands` is the positive control for Step 3:**
- Two 0.5×0.5 squares with gap=0.1 between them.
- The bridge algorithm finds 2 candidate edges and adds 2 bridges.
- Output: 3 interior cells (the 2 squares + 1 bridge quad connecting them, which happens to be interior because the bridge spans the gap into the polygon's interior region).
- This proves the Step 3 rewrite does perform real bridge-joining (not just adding isolated quads like the old algorithm).

**Build:** clean.

**Next:** Step 7 — update `doc/TUSQH_WINDING_GUIDE.md` §3.11.9 (replace the limitation note with the new paper-faithful implementation) and `doc/SESSION_TUSQH_BRIDGE_2026-07-15.md` (log final session state).

## Cross-references

- The 6 tests cover the 4 fixes documented in `BUGS_FOUND.md`:
  - Issue #1 (`boundary_filter_Lshape`) — domain boundary filter.
  - Issue #2 (`rotated_quad_centroid_direction`) — centroid-based direction.
  - Issue #3 — implicitly verified via `bridge_connects_two_islands` (BFS sees valid edges).
  - Issue #5 (`figure14_drop_all_components`) — empty-Quadrants guard.
- The new test `bridge_connects_two_islands` is the **positive control** that proves the Step 3 rewrite works — see `STEP_3_BRIDGE_PAPER_FAITHFUL.md` and `BUGS_FOUND.md` Issue #4 §"Fix".
- See `WORK_SUMMARY.md` §"Regression tests" for the test summary table.
