# Current state — `feature/paper-faithful-bridge`

## Branch

- **Name:** `feature/paper-faithful-bridge`
- **Base:** `develop-felipe` @ commit `d174e6c` ("Union Winding Numbers
  con Templates")
- **Current HEAD:** `d174e6c` (no new commits yet)
- **Working tree:** dirty — all changes are unstaged

## File-by-file diff summary

### Source code

| File | Insertions | Deletions | Purpose |
|------|-----------:|----------:|---------|
| `src/Mesher.cpp` | +1023 | -? | Steps 1–4 + bug fixes |
| `src/Mesher.h` | +152 | -? | Function declarations |
| `src/CMakeLists.txt` | +10 | -? | New sources added |
| `src/Main.cpp` | +57 | -? | New CLI flags (if any) |
| `src/MeshPoint.h` | +74 | -? | API changes |
| `src/MeshPoint.cpp` | +26 | -? | API changes |
| `src/Visualization/VolumeFractionVTKWriter.h` | +24 | -? | Debug output |
| `src/Visualization/VolumeFractionVTKWriter.cpp` | +187 | -? | Debug output |

### Documentation

- **Modified:** `doc/TUSQH_WINDING_GUIDE.md` (+305 lines, §3.11.9
  rewritten).
- **New:** `doc/STEP_0..7_*.md` (8 files), `doc/BUGS_FOUND.md`,
  `doc/WORK_SUMMARY.md`, `doc/FUTURE_WORK.md`,
  `doc/SESSION_TUSQH_BRIDGE_2026-07-15.md`.

### Scripts

- **New:** `scripts/test_bridge_clean_info.py` (6 regression tests).
- **New:** `scripts/analyze_output.py` (quantitative validation).

### Data

- **New:** `data/tusqh_boundary_filter.poly`,
  `data/tusqh_bridge_small_gap.poly`,
  `data/tusqh_rotated.poly`,
  `data/tusqh_figure14.poly`,
  `data/Agua.poly`.

### Build artifacts

- **Untracked:** `build/` directory (binary `mesher_roi`).
- **Untracked:** debug VTKs `_subcell_edge.vtk`, `_subcell_vertex.vtk`,
  `volume_fraction_debug.vtk`, `volume_fraction_debug_samples.vtk`.

## Verification status

| Check | Status | Last verified |
|-------|--------|---------------|
| `cmake --build build` | ✅ Clean | 2026-07-15 |
| Regression tests 6/6 | ✅ Pass | 2026-07-15 |
| Chesapeake Bay quantitative validation | ✅ Pass | 2026-07-15 |
| Visual ParaView verification | ❌ **Not done** — user to verify | — |

## What works

The 7-step plan completed in this branch:

1. ✅ Step 0 — Branch created.
2. ✅ Step 1 — AllOutside cells preserved for BFS.
3. ✅ Step 2 — `isInteriorCell` helper.
4. ✅ Step 3 — Bridge loop rewritten to be paper-faithful.
5. ✅ Step 4 — Output filter drops AllOutside.
6. ✅ Step 5 — Quantitative validation passes.
7. ✅ Step 6 — Regression tests 6/6 pass.
8. ✅ Step 7 — Documentation overhaul complete.

Final Chesapeake Bay metrics:

- 549 interior cells in output
- 110 components
- 125 bridges added in 2 iterations
- 970 AllOutside cells dropped

## What doesn't work yet

These are documented in `FUTURE_WORK.md` and ordered by priority:

1. Visual ParaView verification (P1).
2. Final build + tests re-run after docs are committed (P1).
3. `-B` flag to disable bridges for A/B comparison (P2).
4. `-K0` override for `preRefineForTusqh` base level (P2).
5. 3D pinch templates (P2).
6. BFS performance optimization (P2).
7. Persistence diagrams with Gudhi (P2, requires new dep).
8. Doxygen comments on new functions (P3).
9. CMake flag to gate the new path (P3).
10. `Logger` class replacing scattered `cout` calls (P3).

## Known limitations

The branch makes certain design choices that the user should be
aware of:

1. **Bridge quads are classified via WindingNumberVisitor after
   split.** They start as `Unknown` and need to be re-classified.
   This is documented in `BUGS_FOUND.md` and implemented in
   `Mesher.cpp:3560-3590`.

2. **`Unknown` is NOT considered interior by `isInteriorCell`.**
   Tested and rejected — including it produces infinite loops in the
   BFS. `Unknown` quads are always reclassified before being used in
   the BFS.

3. **BFS is full-grid (interior + exterior preserved), not
   interior-only.** This is correct per paper §3.4. See `BUGS_FOUND.md`
   Issue #4 and `WORK_SUMMARY.md` §3.

4. **`preRefineForTusqh` triggers on `Quadrants.size() <= 2 && segments >= 100`.**
   Hardcoded. Override via P2 in `FUTURE_WORK.md`.

5. **No `-B` flag yet.** Without it, the only way to compare old vs.
   new algorithm is to revert the code.

## Uncommitted changes

The branch has **zero new commits**. All changes are in the working
tree. The user has explicitly not requested a commit.

To commit when ready:

```bash
git add doc/ scripts/ src/ data/  # don't add build/, *.vtk
git status                         # verify
git commit -m "Paper-faithful TUSQH bridge-joining with sub-cell VF"
```

## Review checklist

Before merging to `develop-felipe`, verify:

- [ ] Visual ParaView check (user).
- [ ] Build is clean from fresh checkout.
- [ ] 6/6 tests still pass.
- [ ] Chesapeake Bay metrics are unchanged from
      `WORK_SUMMARY.md` table.
- [ ] No new bugs introduced in the documentation commits.
- [ ] All untracked files (`build/`, debug VTKs) are properly
      gitignored.

## See also

- `WORK_SUMMARY.md` — high-level summary of what was done.
- `BUGS_FOUND.md` — every bug and how it was fixed.
- `FUTURE_WORK.md` — outstanding work, prioritized.
- `doc/STEP_0..7_*.md` — detailed per-step reports.
- `doc/TUSQH_WINDING_GUIDE.md` §3.11.9 — user guide, rewritten.
- `doc/SESSION_TUSQH_BRIDGE_2026-07-15.md` — session log.
