# Step 2 — Add `isInteriorCell` helper

**Status:** ✅ SUCCESS

**Files modified:**
- `src/Mesher.h` — added static method declaration after `isEdgeOnDomainBoundary`
- `src/Mesher.cpp` — added inline definition after `computeExteriorDirection`

**Declaration:**
```cpp
// Paper-faithful helper: returns true iff a Quadrant is an
// "interior" cell of the cubical complex, i.e. it has positive
// volume fraction (AllInside) or is on the boundary (Mixed).
static bool isInteriorCell(const Quadrant &q);
```

**Definition:**
```cpp
bool Mesher::isInteriorCell(const Quadrant &q) {
    switch (q.getWindingState()) {
        case WindingState::AllInside:
        case WindingState::Mixed:
            return true;
        case WindingState::AllOutside:
        case WindingState::Unknown:
        default:
            return false;
    }
}
```

**Rationale:** Centralises the "is this a real mesh cell?" check used by the bridge-joining BFS, `saveOutputMesh`, `linkElementsToNodes`, `detectFeatureQuadrants`, `detectInsideNodes`, etc. Without this helper, each consumer would have to spell out the `WindingState` switch.

**Build:** clean.

**Verification:** none yet — the helper is unused until Step 3 wires it into the bridge-joining loop. Sanity check: builds without warnings.

**Next:** Step 3 — rewrite the bridge-joining loop in `resolveArchipelagos` to use the new helper and find shared edges between interior components.

## Cross-references

- `isInteriorCell` deliberately returns **false** for `WindingState::Unknown`. Tested and rejected: including `Unknown` in the BFS produces infinite loops (revisit during bridge reclassification). See `BUGS_FOUND.md` Issue #4 §"Fix" for the full reasoning.
- The helper is used in `STEP_3_BRIDGE_PAPER_FAITHFUL.md` (BFS filter) and `STEP_4_OUTPUT_FILTER.md` (output filter).
- See `WORK_SUMMARY.md` §"Design choices" for the broader design pattern.
