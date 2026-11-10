## What this pull request does

Scenario or component ID: <!-- e.g. T1.4 or E2 -->
Linked issue(s): <!-- claim issue, result issue -->

## Checklist (TEAM.md, section 11)

- [ ] One component or one scenario only (about 400 changed lines at most, generated output excepted)
- [ ] Pre-registration committed in `team/notes/<ID>.md` before production runs (scenarios only)
- [ ] T0 passes on this engine commit (CI)
- [ ] Results written by scripts under `results/<ID>/<run-id>/` with `stamp.json` and `summary.json`
- [ ] Conservation drift and convergence (dt/2, dt/4, particle counts) reported
- [ ] Every statement labelled: mathematics / measurement (with source) / finding of the virtual universe (law Lx)
- [ ] Result compared with the prediction
- [ ] No edits to verified formulas, verification scripts, sources or pass criteria (unless that is the task and the lead approved)
- [ ] Code changed? `ARCHITECTURE.md` updated: registry, Mermaid diagrams, contracts, change log (D18); `python tools/check_architecture.py --base origin/main` passes
- [ ] `team/LOG.md`, `team/FINDINGS.md` and `docs/OPEN-SCENARIOS.md` updated

## For AI agents

Agent / model / version:
Human sponsor:
What I could not verify:
