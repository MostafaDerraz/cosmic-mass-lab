# Virtual universe of cosmic masses: the web platform

`simulator/web/` is the interactive universe (tasks UI1-UI4 of `docs/WORK-PLAN.md`). It runs the **unchanged, verified engine** `cosmicmass` inside the browser with [Pyodide](https://pyodide.org) (Python 3.12, NumPy 1.26.4, the pinned versions). There is no second engine in JavaScript: every number on the page comes from `cosmicmass.runner.run`.

## Online

Published with GitHub Pages from `main`: <https://mostafaderraz.github.io/cosmic-mass-lab/simulator/web/>. The empty file `.nojekyll` at the repository root is needed: without it, Pages hides every file whose name starts with `_`, such as the engine's `__init__.py`.

## Run it locally

From the repository root:

```bash
python -m http.server 8765
```

Then open <http://localhost:8765/simulator/web/>. The first load takes about 20 seconds (Python and the engine are downloaded into the browser); a two-body scene runs in 10 to 40 seconds.

## What it does (UI1, first version)

- Scenes: hydrogen-like pair (demo T1.4), Earth and Moon, two protons from rest, Earth and Moon falling from rest; or a custom universe of up to 16 bodies.
- Laws: L0 (reference) and L1a (hypothesis) with its strength, sign and terms; the frame angle; an overlay of the same start under L0.
- Animated trajectories in the Real plane, with a time slider; readouts of the steps, the energy change at the end **and** its largest excursion (D21), momentum, mass and charge drifts.
- Download of the `settings.json` (runs the same with `python -m cosmicmass run settings.json`) and of the result.

## Files

| File | Role |
|---|---|
| `index.html`, `style.css`, `app.js` | the page: editors, canvas, readouts |
| `worker.js` | loads Pyodide and the engine files in a background thread, runs scenes |
| `bridge.py` | Python side: calls the engine, returns the trajectory as JSON |
| `engine-files.json` | the list of engine files the browser loads; `engine/tests/test_web_universe.py` checks that it matches the engine |

When an engine module is added, regenerate the list:

```bash
python -c "import json,pathlib; fs=sorted(p.as_posix().split('engine/',1)[1] for p in pathlib.Path('engine/cosmicmass').rglob('*.py')); json.dump({'root':'../../engine/','files':fs}, open('simulator/web/engine-files.json','w'), indent=1)"
```

## Not yet (next UI tasks)

- The Composer plane (each body's `mu` drawn in the leaning frame) and live readings per body (UI1, second part).
- Compare mode with several laws side by side (UI3); shareable links and saved runs (UI4); publication through GitHub Pages.
- Speed: Pyodide is several times slower than native Python; scenes with many bodies need the fast solver (E2.5).
