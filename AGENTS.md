# Instructions for AI agents

AI agents (coding assistants, autonomous agents, research agents) are welcome contributors. Agents follow every rule in `CONTRIBUTING.md`, plus the rules below.

## Before you start

1. Read these files:
   - `CONTRIBUTING.md`;
   - `docs/03-report-verified-formulas.md`;
   - `docs/formulas-compendium/formulas-compendium.pdf` (or `.tex`);
   - `docs/05-scenarios-and-tests.md`.
2. Run every script in `verification/` (commands in `.github/workflows/verify.yml`). If any fails, stop and report it; do not work around it.
3. Work on one scenario at a time, and only one that is claimed for you on the scenario board.

## While you work

4. **Use only the verified formulas.** Do not derive new physics silently. If a scenario needs a formula that is not in the compendium, write it in your pull request as a **derivation**, with a symbolic check (SymPy) that it follows from the compendium.
5. **Never edit:**
   - files in `verification/`, unless your task is a verification script;
   - files in `sources/`;
   - the pass criteria of a test.
6. **Never present a hypothesis as established.** Label every L1-L8 result "finding of the virtual universe (law Lx)".
7. **Never invent numbers, citations or measurements.** Every measured value needs a source that a human can check.
8. **Show your work.** Commit the script, the settings stamp, the raw output and the convergence check (T7). Results that exist only in prose are not accepted.
9. **Be conservative with compute.** Start with the smallest setup that answers the question, and scale up only after it passes.

## When you submit

10. **Pull request title:** `[agent] T<id>: <short description>`.
11. **Pull request body:**
    - agent name, model and version;
    - the human sponsor (GitHub handle);
    - the test ID and its status: pass, fail or finding;
    - the T0 status of the engine commit used;
    - what you could not verify.
12. **A second, independent reproduction is required** (another agent or a human) before a finding is marked "confirmed".
13. **If you are unsure, say so.** An honest "unclear" is better than a confident error.
