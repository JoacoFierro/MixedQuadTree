# Documentation index

This directory contains the documentation for the
`feature/paper-faithful-bridge` branch. Start here.

## Quick links

| If you want to… | Read this |
|-----------------|-----------|
| Know the goal and final results | [`WORK_SUMMARY.md`](WORK_SUMMARY.md) |
| Understand a specific bug and how it was fixed | [`BUGS_FOUND.md`](BUGS_FOUND.md) |
| See what's still pending | [`FUTURE_WORK.md`](FUTURE_WORK.md) |
| Check the branch state before merging | [`CURRENT_STATE.md`](CURRENT_STATE.md) |
| Read about a single implementation step | [`STEP_0..7_*.md`](STEP_0_BRANCH_SETUP.md) |
| Use the API (`-T -J -K -F -L` flags, internal functions) | [`TUSQH_WINDING_GUIDE.md`](TUSQH_WINDING_GUIDE.md) |
| Read the original session log | [`SESSION_TUSQH_BRIDGE_2026-07-15.md`](SESSION_TUSQH_BRIDGE_2026-07-15.md) |

## Document structure

```
doc/
├── README.md                       (this file — index)
│
│   ── Master documents (read first) ──────────────────────────
├── WORK_SUMMARY.md                 goal + plan + final metrics
├── BUGS_FOUND.md                   every bug with fix
├── FUTURE_WORK.md                  TODOs P1-P5
├── CURRENT_STATE.md                branch status, what's ready
│
│   ── Implementation reports (read for details) ──────────────
├── STEP_0_BRANCH_SETUP.md          branch creation
├── STEP_1_PRESERVE_EXTERIOR.md     preserve AllOutside cells
├── STEP_2_IS_INTERIOR_CELL.md      helper function
├── STEP_3_BRIDGE_PAPER_FAITHFUL.md rewrite of bridge loop
├── STEP_4_OUTPUT_FILTER.md         filter exterior from output
├── STEP_5_VALIDATION.md            Chesapeake Bay metrics
├── STEP_6_REGRESSION_TESTS.md      6/6 tests passing
├── STEP_7_DOCUMENTATION.md         this documentation overhaul
│
│   ── Reference documents (use as needed) ────────────────────
├── TUSQH_WINDING_GUIDE.md          API guide, flags, call-graph
├── SESSION_TUSQH_BRIDGE_2026-07-15.md  original session log
│
└── ... (other legacy docs if any)
```

## Reading order

For a new contributor or reviewer joining the branch, we recommend
reading in this order:

1. **[`WORK_SUMMARY.md`](WORK_SUMMARY.md)** — 5 min read. Get the
   high-level picture: what was built, why, what worked, what's
   pending.
2. **[`CURRENT_STATE.md`](CURRENT_STATE.md)** — 2 min read. Know the
   branch state: uncommitted, ready for review, 6/6 tests pass.
3. **[`BUGS_FOUND.md`](BUGS_FOUND.md)** — 10 min read. Skim the 7
   issues to understand what was non-trivial.
4. **[`STEP_3_BRIDGE_PAPER_FAITHFUL.md`](STEP_3_BRIDGE_PAPER_FAITHFUL.md)**
   — 10 min read. The heart of the implementation: the new bridge
   loop.
5. **[`TUSQH_WINDING_GUIDE.md`](TUSQH_WINDING_GUIDE.md)** §3.11.9 —
   5 min read. The user-facing description of the new pipeline.
6. **[`FUTURE_WORK.md`](FUTURE_WORK.md)** — 5 min read. Know what's
   pending so you can pick something up.

For someone picking up a specific task (e.g., adding the `-B` flag):

1. **[`FUTURE_WORK.md`](FUTURE_WORK.md)** §2.1 — the task description.
2. **[`BUGS_FOUND.md`](BUGS_FOUND.md)** Issue #4 — context on the
   bridge-joining rewrite that the `-B` flag would disable.
3. **[`STEP_3_BRIDGE_PAPER_FAITHFUL.md`](STEP_3_BRIDGE_PAPER_FAITHFUL.md)**
   — the code location (`Mesher.cpp:3361-3520`).
4. **[`TUSQH_WINDING_GUIDE.md`](TUSQH_WINDING_GUIDE.md)** §4 — the
   CLI flag parser location (`src/Main.cpp`).

## Cross-references

Most documents cross-link to others. The most-cited are:

- `WORK_SUMMARY.md` — cited from `BUGS_FOUND.md`, `FUTURE_WORK.md`,
  `CURRENT_STATE.md`, all `STEP_*.md`, `TUSQH_WINDING_GUIDE.md` §0,
  `SESSION_TUSQH_BRIDGE_2026-07-15.md`.
- `BUGS_FOUND.md` — cited from `STEP_3`, `STEP_4`, `STEP_6`,
  `CURRENT_STATE.md`.
- `FUTURE_WORK.md` — cited from `CURRENT_STATE.md`, all `STEP_*.md`.
- `CURRENT_STATE.md` — cited from `STEP_7`, `WORK_SUMMARY.md`.

## Conventions

- **File names** in `doc/` are uppercase with underscores
  (`STEP_3_BRIDGE_PAPER_FAITHFUL.md`).
- **Code references** are `file:line` (e.g., `Mesher.cpp:3361`).
- **Issue references** are `#N` (e.g., `Issue #4`).
- **Step references** are `Step N` (e.g., `Step 3`).
- **Status markers** are ✅ (done), ⚠️ (partial / known issue),
  ❌ (not done / not applicable).

## Generated vs hand-written

All files in `doc/` are hand-written Markdown. There are no
auto-generated docs in this repo (Doxygen comments are sparse; see
`FUTURE_WORK.md` §3.1).

## Last updated

2026-07-15.

For questions, contact the branch author or file an issue.
