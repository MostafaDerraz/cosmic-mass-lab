# Decisions

Numbered, binding decisions. To propose one, add it with status "proposed" and open an issue. The maintainer approves or rejects it. Never delete a decision; supersede it with a new one.

| ID | Decision | Status | Date |
|---|---|---|---|
| D1 | The paper's law L0 (`Re(mu) dv/dt = Re[mu(E_cal + v x B_cal)]` with F1-F4) is the reference law. Every other reading (L1-L8) is a hypothesis. | approved | 2026-11-10 |
| D2 | Engine, outputs and documents use only the language of the cosmic mass and General Numbers. | approved | 2026-11-10 |
| D3 | T0 is a gate: no coupling result counts unless the engine commit used passes T0. | approved | 2026-11-10 |
| D4 | T9 (the measured limits of the book's Chapter 10) must pass before exploring L1-L8. | approved | 2026-11-10 |
| D5 | Verified formulas (`docs/formulas-compendium/`) and verification scripts are not edited to make results pass. | approved | 2026-11-10 |
| D6 | The book's Chapters 7 (Standard Model) and 8 (Closure Law) are not used as engine inputs; their physical claims are unverified. | approved | 2026-11-10 |
| D7 | Units: SI is the reference; scaled units (`G_N = K = 1`, so `s = 1`) are allowed for visibility, with the SI run shown or the reason it cannot be run. | approved | 2026-11-10 |
| D8 | Licence: MIT for code, CC BY 4.0 for documents and findings; the sources remain the author's work. | approved | 2026-11-10 |
| D9 | The repository is public. | approved | 2026-11-10 |
| D10 | Default numerics: velocity Verlet (Boris with `B_cal`), adaptive time step, softening 1e-3 of the smallest initial separation. | approved | 2026-11-10 |
