# Step 1 — Preserve `AllOutside` cells in `windingSubdivide`

**Status:** ✅ SUCCESS (with note about Step 4 dependency)

**File modified:** `src/Mesher.cpp` lines ~2030-2068

**Change:** The TUSQH post-processing loop no longer discards `AllOutside` cells. Instead, it preserves them in `Quadrants` with their `WindingState` intact, so the cubical complex remains intact for the bridge-joining algorithm.

```cpp
case WindingState::AllOutside:
    q.getIntersectedEdges().clear();
    Quadrants.push_back(std::move(q));
    ++preservedOutside;  // was: ++droppedOutside; break;
    break;
```

**Rationale:** The paper's bridge-joining algorithm (§3.4) iterates the background grid to find edges between two interior components. If we drop AllOutside cells, MapEdges loses those edges and the topology is broken. Preserving them lets the BFS see all interior-adjacent pairs.

**Verification (Chesapeake Bay, `-a 3 -T -J -K 2 -F 0.5 -L 1`):**
- Before: `TUSQH dropped 345 AllOutside cells`
- After: `TUSQH preserved 345 AllOutside cells in cubical complex`
- Bridge-joining: still 0 bridges added (expected — the new bridge algorithm is Step 3).
- 74 components (down from 81 with original algorithm, same as pre-fix).

**Issue observed:** `saveOutputMesh` does NOT filter exterior cells, so the current VTK output now contains ALL cells including `AllOutside`. Output count went from 645 (pre-fix) to 1019 (post Step 1). This is fixed in Step 4 by adding an `isInteriorCell` filter inside `saveOutputMesh`.

**Regression tests:** NOT yet run — they will be updated in Step 6 after the full pipeline is in place.

**Next:** Step 2 — add `isInteriorCell` helper to `Mesher.h`.

## Cross-references

- The `AllOutside` preservation is a prerequisite for Step 3's full-grid BFS — without it, `MapEdges` would be missing the edges between exterior cells, breaking the topology. See `STEP_3_BRIDGE_PAPER_FAITHFUL.md` for the algorithm that consumes the preserved cells.
- This step fixes part of Issue #4 in `BUGS_FOUND.md` (the highest-severity bug, "bridge-joining produces isolated quads").
- See `WORK_SUMMARY.md` §3 for the rationale of "preserve-then-filter" pattern.
