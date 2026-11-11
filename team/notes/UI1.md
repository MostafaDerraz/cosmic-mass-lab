# UI1: the web universe, first version

- **Owner:** Claude Opus 5.5 in Claude Code, local (sponsor @MostafaDerraz)
- **Claimed:** 2026-11-11, asked by the maintainer ("start a web platform with a frontend of this")
- **Depends on:** version 0.1.0 (merged)
- **Branch:** UI1-web-universe

## Goal
A web page where anyone can build a small universe of cosmic masses, run it under L0 or L1a, and watch it, with every number computed by the verified engine.

## Design (UI0, decision D22 proposed)
- The browser runs the **unchanged engine** with Pyodide 0.26.4 (Python 3.12, NumPy 1.26.4 = the pin, mpmath, SymPy) in a web worker, so the page stays responsive.
- `bridge.py` calls `cosmicmass.runner.run` with `out_root=None` and returns at most 1500 frames. In the browser there is no git, so the stamp's engine commit is `web-v<version>`; nothing else is changed.
- Alternatives rejected: a JavaScript copy of the engine (a second engine that could disagree with the verified one); a server (costs, hosting, accounts).

## Progress (newest first)
- 2026-11-11: first version works in the browser: engine ready in about 20 s; hydrogen-like pair under L1a with the L0 overlay in 33 s (10,034 steps), the ellipse and the circle of T1.4. `engine/tests/test_web_universe.py`: the manifest equals the engine's files, and the bridge runs an Earth-moon arc natively (energy change below 1e-9).

## Results
Software only; no finding. The demo picture reproduces T1.4 (exploratory).

## Open questions
- Publication with GitHub Pages: a repository setting for the maintainer.
- Speed: the adaptive step limits every body to 1e-3 of its nearest-neighbour distance using its **absolute** speed, so scenes that move as a whole (for example Sun, Earth and Moon) take millions of steps. Worth a scenario (relative speed instead), since it is a numerical convention.
- Next: the Composer plane, per-body readings, shareable links (UI1 second part, UI3, UI4).

## Handoff
- State: first version on `UI1-web-universe`, pull request to open.
- Exact next step: review; enable GitHub Pages; then the Composer plane.
