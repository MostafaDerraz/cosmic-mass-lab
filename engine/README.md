# Engine

Components E1-E9: tasks in `docs/WORK-PLAN.md`, modules and interface contracts in `ARCHITECTURE.md`. The package skeleton exists (P1): `pip install -e .[dev]` from the repository root, then `pytest`. Implemented so far: the algebra core E1 (`constants.py`, `frame.py`, `gn.py`, `gnarray.py`, `readings.py`, `matter.py`). Rule D18: every code change updates `ARCHITECTURE.md` in the same pull request.
