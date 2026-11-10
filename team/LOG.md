# Project log

Newest entries first. One entry per work session, even if nothing worked. Format:

```
### YYYY-MM-DD - @handle (or agent / model, sponsor @handle) - <scenario ID or topic>
- Did:
- Result:
- Next:
- Blocked by / questions:
```

---

### 2026-11-10 - Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2) (sponsor @MostafaDerraz) - MVP-4 (E5.1, E5.3 minimal, E4.1, E4.2)
- Did: wrote `engine/cosmicmass/fields/sample.py` (new module: `FieldSample`, `FieldSolver`), `state.py`, `laws/base.py`, `laws/l0.py` and a minimal `settings.py`, with `test_state_laws.py` and `test_settings.py` (51 new tests). Opened the notebook `team/notes/MVP-4.md`; updated the board and the MVP checklist. (MVP-3, the CI job, is on its own branch `P2-ci-tests`, with its own log entry there.)
- Result: `pytest` gives 115 passed, coverage 100%. Two bodies under L0 feel `(G_N m m' - K q q')/r^2` to 1e-15, with equal and opposite forces; the L0 force is bit-for-bit the same at every frame angle. Mathematics from the compendium; no finding.
- Next: MVP-5 (`fields/direct.py`: static direct sum, softening, backgrounds) on `E2-fields`, cut from this branch.
- Blocked by / questions: none. Choices to review are in the notebook (`Law` as a base class; `LawSpec.moving` as a string; the new module `fields/sample.py`).
- Architecture: updated (sections 2, 4, 5, 10, 11)

### 2026-11-10 - Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2) (sponsor @MostafaDerraz) - MVP-2 (E1.3 to E1.6)
- Did: wrote `engine/cosmicmass/gnarray.py` (E1.3), `readings.py` (E1.4), `matter.py` (E1.5) and the tests `test_gnarray.py`, `test_readings_matter.py`, `test_compendium.py` (E1.6); 36 new tests. Updated the notebook `team/notes/E1.md`, the board (E1 `done`) and the MVP checklist.
- Result: `pytest` gives 64 passed, coverage 100% of `engine/cosmicmass`. Arrays agree with scalars to 1e-15 on 10,000 values; the product `mu E_cal` matches the compendium in both coordinates; `Re(mu_a mu_b)` is the same number at every frame angle. Mathematics from the compendium; no finding of the virtual universe.
- Next: MVP-3 (P2, CI runs pytest) on `P2-ci-tests`, cut from `main`; then MVP-4 on `E4-l0`, cut from this branch.
- Blocked by / questions: none new. `Com(z)` returns the coordinate `b`, and `CosmicMass.mu(frame)` takes the frame (see the notebook).
- Architecture: updated (sections 2, 4, 5, 10, 11)

### 2026-11-10 - Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2) (sponsor @MostafaDerraz) - MVP-1 (E1.1, E1.2)
- Did:
  - read `docs/MVP-PLAN.md` and decision D20; took the MVP sequence as one claim;
  - wrote `engine/cosmicmass/constants.py` and `frame.py` (E1.1) and `gn.py` (E1.2), with `engine/tests/test_constants_frame.py` and `test_gn.py` (25 new tests);
  - opened the notebook `team/notes/E1.md`; updated the board, the MVP checklist and `engine/README.md`.
- Result: `pytest` gives 28 passed, coverage 100% of `engine/cosmicmass`. Rule of the Square exact at 50 digits; `z * (1/z) = 1`; products and quotients agree with the bridge to C. All of it is mathematics from the compendium; no finding of the virtual universe.
- Next: MVP-2 (`gnarray.py`, `readings.py`, `matter.py`, `test_compendium.py`) on `E1b-arrays-readings`, cut from this branch.
- Blocked by / questions: three open questions in `team/notes/E1.md`: the value of `K` differs between `docs/05` (CODATA 2018) and the verified constants (CODATA 2022, used here); `c` and `e` in scaled units are not fixed by D7 (assumed Stoney units, `c = e = 1`); the compendium has no equation numbers (tests cite paragraph titles).
- Architecture: updated (sections 2, 4, 5, 10, 11)

### 2026-11-10 - Claude Opus 5.5 in Claude Code, local (sponsor @MostafaDerraz) - P1 review and delivery
- Did: applied the P1 patch produced by the claude.ai chat session below (it had no push access), reviewed it as a different model (D14), and re-ran the done-when check on Windows: clean clone of the branch, fresh Python 3.12.10 environment, `pip install -e .[dev]` then `pytest` gives 3 passed (100% coverage); `python tools/check_architecture.py --base origin/main` passes. Corrected the session label in this log and in `team/notes/P1.md`. Pushed the branch and opened the pull request.
- Result: P1 works on Linux (producer) and Windows (reviewer). Remarks for the reviewer of record: the dev extras `pytest` and `pytest-cov` are not pinned (resolved to 9.2.0 and 7.1.0); suggest pinning them in P2 together with a check that `pyproject.toml` and `requirements.txt` pins agree.
- Next: maintainer review and merge (D17); then P2.
- Blocked by / questions: none.
- Architecture: no code changed in this session (the architecture update is in the P1 commit)

### 2026-11-10 - Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2) (sponsor @MostafaDerraz) - P1
- Did:
  - ran the six verification scripts in the pinned environment (Python 3.12.3): 106, 204, 253 checks, 144 laws, 9,304 Gazette checks, 79 General Sinus checks, all pass; `python tools/check_architecture.py` passed before starting;
  - added `pyproject.toml` (package `cosmicmass` from `engine/cosmicmass/`, Python >= 3.12, the three pins of `requirements.txt`, dev extras `pytest` and `pytest-cov`, pytest pointed at `engine/tests`);
  - added `engine/cosmicmass/__init__.py` (`__version__ = "0.0.1"`) and `engine/tests/test_import.py` (3 tests);
  - added packaging and test leftovers to `.gitignore`; updated `engine/README.md`; opened the notebook `team/notes/P1.md`; marked P1 `done` on the board.
- Result: on a clean copy of the branch in a fresh Python 3.12 environment, `pip install -e .[dev]` then `pytest` gives 3 passed. No physics code and no finding; `team/FINDINGS.md` is unchanged.
- Next: P2 (CI job for tests and coverage) can start; it should also keep `pyproject.toml` and `requirements.txt` pins in step. P4 and P7 are unblocked too.
- Blocked by / questions: no claim issue was opened for P1 (the task was assigned directly by the sponsor, who is the maintainer). The pull request needs a review by someone other than the author (D17).
- Architecture: updated (sections 2, 4, 11)

### 2026-11-10 - @MostafaDerraz with Claude Opus 5.5 (AI assistant) - project setup
- Did:
  - verified the three papers and the book's algebra (106, 204, 253 checks, 144 laws, 9,304 Gazette checks, 79 General Sinus checks);
  - audited General Sinus and corrected five errors (see `docs/03-report-verified-formulas.md`);
  - wrote the formulas compendium, the test plan T0-T9, the engine architecture and the general note;
  - created this repository with CI verification, opened it to the public, and wrote the rules (`CONTRIBUTING.md`, `AGENTS.md`), the team guide (`TEAM.md`) and the open scenarios board.
- Result: the foundation is verified and CI is green. No simulator code exists yet.
- Next: milestone 1. Engine components E1-E4, then the T0 gate.
- Blocked by / questions: none. Engine lead wanted.
