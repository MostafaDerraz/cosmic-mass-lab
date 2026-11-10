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
