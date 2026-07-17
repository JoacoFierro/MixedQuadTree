# Step 4 — Output filtering (interior cells only)

**Status:** ✅ SUCCESS

**Files modified:** `src/Mesher.cpp` (drop logic in `resolveArchipelagos` lines ~3646-3690)

**Changes:**

The post-archipelago output (`_postarchipelago.vtk`) is now correctly filtered to **interior cells only**. Previously the algorithm preserved AllOutside cells in `Quadrants` for the BFS / bridge-joining topology (Step 1), but the post-resolve save-output step would save all of them too — yielding a mesh that included zero-VF cells around every island.

**Filter applied in `keepQuad[]`:**
```cpp
if (compOfQuad[qi] >= 0 &&
    compSize[compOfQuad[qi]] >= (int)minComponentCells &&
    isInteriorCell(Quadrants[qi])) {
    keepQuad[qi] = true;
}
```

The third condition (`isInteriorCell`) drops AllOutside cells from the output regardless of which component they belong to. They were kept in `Quadrants` only so the BFS could route through them and the bridge algorithm could identify edge pairs that straddle exterior corridors.

**Verification (Chesapeake Bay, `-a 3 -T -J -K 2 -F 0.5 -L 1`):**

| Metric | Before Step 4 | After Step 4 |
|--------|---------------|--------------|
| `_tusqh.vtk` cells | 1019 | 1019 (unchanged — debug output, both interior + exterior) |
| `_postarchipelago.vtk` cells | 1519 (incl. 970 exterior) | 549 (interior only) |
| `_postarchipelago.vtk` AllOutside cells | 970 | 0 |

The diagnostic print now correctly shows `970 small-component quads dropped` (these are the AllOutside cells that were preserved for the BFS but dropped from the output).

**Important:** The `_tusqh.vtk` debug file (saved BEFORE `resolveArchipelagos`) still contains both interior and exterior cells — this is intentional and useful for debugging the TUSQH winding classification.

**Build:** clean.

**Next:** Step 5 — visual validation (ParaView screenshots / analysis of the post-archipelago output to confirm the algorithm is doing what we expect for Chesapeake Bay).

## Cross-references

- The output filter uses the `isInteriorCell` helper from `STEP_2_IS_INTERIOR_CELL.md` and is the consumer of the `AllOutside` cells preserved in `STEP_1_PRESERVE_EXTERIOR.md`. Without either, this step would not be possible.
- The filter is the last piece of the "preserve-then-filter" pattern documented in `WORK_SUMMARY.md` §3.
- See `BUGS_FOUND.md` Issue #5 for the `saveOutputMesh` empty-Quadrants guard, which is a related crash protection.
- See `STEP_5_VALIDATION.md` for the quantitative analysis confirming this filter correctly drops 970 `AllOutside` cells from the Chesapeake Bay output.
