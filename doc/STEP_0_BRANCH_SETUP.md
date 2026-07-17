# Step 0 — Branch setup

**Status:** ✅ SUCCESS

**Date:** 2026-07-15

**Action:**
```bash
git checkout -b feature/paper-faithful-bridge
```

**Starting point:** `d174e6c Union Winding Numbers con Templates` on `develop-felipe`.

**Branch:** `feature/paper-faithful-bridge` (new local branch, not pushed).

**Pre-existing changes carried over** (uncommitted on `develop-felipe`):
- `src/Mesher.cpp`, `src/Mesher.h`, `src/Main.cpp`, `src/MeshPoint.{h,cpp}`,
  `src/Visualization/VolumeFractionVTKWriter.{h,cpp}`, `src/CMakeLists.txt`
- New: `src/SubgridSampler.{h,cpp}`, `src/SubcellVFData.h`,
  `scripts/test_bridge_clean_info.py`,
  `data/{tusqh_boundary_filter,tusqh_figure14,tusqh_rotated,tusqh_bridge_small_gap}.poly`,
  `data/Agua.poly`, `doc/SESSION_TUSQH_BRIDGE_2026-07-15.md`

These are the "Pre-fix #4" state — the starting point for this work.

**Next:** Step 1 — preserve `AllOutside` cells in `windingSubdivide`.

## Cross-references

- See `WORK_SUMMARY.md` for the high-level plan and final results.
- See `BUGS_FOUND.md` for the bugs fixed during this work.
- See `FUTURE_WORK.md` for outstanding work.
- See `CURRENT_STATE.md` for the branch's current state (uncommitted, working tree only).
- See `SESSION_TUSQH_BRIDGE_2026-07-15.md` for the original session log that motivated this branch.
