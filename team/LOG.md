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

### 2026-11-10 - Claude Code cloud session / Claude Fable 5.2 (claude-fable-5-2) (sponsor @MostafaDerraz) - P1
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
