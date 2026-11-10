# How the team works

This is the operating guide for everyone on the project: human contributors and AI agents. Read it fully before your first contribution. The rules are in `CONTRIBUTING.md` (everyone) and `AGENTS.md` (agents). This guide explains how we work together day to day.

---

## 1. The objective

**Build a virtual universe made only of cosmic masses (`mu = m + q G`), and use it to find out what the unread second coordinate of their product would do in nature, under each possible reading.**

Concretely, in order:

1. **An engine that passes the consistency gate T0.** Under the paper's law L0 it must reproduce Newton plus Coulomb exactly, conserve mass, charge, energy and momentum, and give results that do not depend on the frame angle.
2. **Reproduce what is already known (T9).** That means the measured limits of the book's Chapter 10, for example a neutron effective charge of 9.0e-19 e under the radial reading.
3. **Explore the readings L1-L8 (T1-T8).** Record for each one what it does: conserving or not, angle-dependent or not, compatible with measurement or not.
4. **Publish the findings**, with credit to everyone who produced and confirmed them.
5. **Deliver an interactive universe** that anyone can use to repeat the experiments.

**What success looks like:** a clear, reproducible table that says, for each reading of the second coordinate, what kind of universe it gives and whether our own universe allows it. Negative answers count as success.

---

## 2. Roles

| Role | Who | Does |
|---|---|---|
| **Maintainer** | Mostafa Derraz | sets the direction, approves decisions (`team/DECISIONS.md`), merges pull requests, owns the theory |
| **Engine lead** | first volunteer on milestone 1 | keeps the engine coherent; reviews every engine change |
| **Scenario owner** | whoever claims a scenario | runs it from claim to result, keeps its notes up to date |
| **Reviewer** | any contributor who did not produce the result | reproduces results independently; marks findings "confirmed" |
| **AI agent** | any agent, with a human sponsor | works like a scenario owner or reviewer under `AGENTS.md` |
| **Scribe** (rotating, weekly) | any contributor | writes the weekly summary in `team/LOG.md` |

One person can hold several roles, but never **producer and reviewer of the same result**.

---

## 3. The shared memory: where everything is written

Work that is not written down did not happen. The team's memory lives in the repository:

| File or folder | What goes there | Who writes |
|---|---|---|
| `team/LOG.md` | the project journal: one short entry per work session, newest first | everyone, every session |
| `team/DECISIONS.md` | every decision that affects others (conventions, laws, units, scope), numbered D1, D2, ... | proposed by anyone, approved by the maintainer |
| `team/FINDINGS.md` | the index of results: scenario, law, result, status (done or confirmed), link | scenario owners and reviewers |
| `team/notes/<ID>.md` | the working notebook of one scenario or component (from `team/notes/_TEMPLATE.md`) | its owner, continuously |
| `docs/OPEN-SCENARIOS.md` | the status board: open, claimed, done, confirmed, blocked | updated with every claim and result |
| issues | claims, results, proposals, questions | everyone |
| pull requests | code, scripts, outputs, notes | everyone |

---

## 4. Before you start: the reading order

Read in this order. Agents must read all of it; humans should.

1. `README.md`: what is in the repository.
2. `TEAM.md`: this guide.
3. `CONTRIBUTING.md`, and `AGENTS.md` if you are an agent.
4. `team/LOG.md`: the last ten entries, so you know what happened recently.
5. `team/DECISIONS.md`: all of it. Decisions are binding.
6. `team/FINDINGS.md`: what is already known.
7. `docs/03-report-verified-formulas.md` and `docs/formulas-compendium/`: the verified foundation.
8. `docs/05-scenarios-and-tests.md`: the full test plan.
9. `team/notes/<ID>.md` of the scenario you want to take, if it exists, so you continue from where someone stopped.

Then run the verification scripts (see `.github/workflows/verify.yml`). If any fails, stop and open an issue.

---

## 5. The work cycle

Every piece of work follows the same seven steps.

### Step 1: Choose

- Pick an `open` scenario from `docs/OPEN-SCENARIOS.md` whose dependencies are `done` or `confirmed`.
- Prefer the lowest milestone that still has open work.

### Step 2: Claim

- Open a "Claim a scenario" issue.
- Mark the scenario `claimed (@handle)` on the board, in a small pull request or by asking the maintainer.
- A claim expires after **14 days without a log entry**.

### Step 3: Open the notebook

- Copy `team/notes/_TEMPLATE.md` to `team/notes/<ID>.md`, or continue the existing one.
- Write your plan before you write code.

### Step 4: Work, writing as you go

- Work on a branch named `<ID>-<short-name>`, for example `T1.4-hydrogen-pair`.
- Update the notebook after every session: what you tried, what happened, what you will do next.
- If you need a convention or choice that affects others, **do not decide alone**. Propose it in `team/DECISIONS.md` (status "proposed") and continue with a clearly marked assumption.

### Step 5: Leave a log entry (every session, even if nothing worked)

Add an entry at the top of `team/LOG.md`:

```
### 2026-11-12 - @handle (or agent name / model, sponsor @handle) - T1.4
- Did: ...
- Result: ...
- Next: ...
- Blocked by / questions: ...
```

### Step 6: Deliver

Open a pull request containing:

- the script;
- the settings stamp;
- the raw output;
- the convergence check;
- the updated notebook;
- a "Scenario result" issue linked to it.

Add a line to `team/FINDINGS.md` with status `done`, and mark the board `done`.

### Step 7: Hand off or close

- **If you stop before finishing:** fill the **Handoff** section of the notebook (state, what works, what does not, the exact next step), write a log entry, and release the claim on the board (`open`, with "see notes").
- **When a reviewer reproduces your result:** they mark it `confirmed` in `team/FINDINGS.md` and on the board, and add their name.

**Definition of done** for a scenario:

- the script runs from a clean checkout;
- T0 passes on the engine commit used;
- convergence is shown;
- the notebook, log, findings index and board are all updated;
- the pull request is merged.

---

## 6. Writing well

- **Label every statement** as mathematics (verified), measurement (cite the source) or finding of the virtual universe (law Lx). Never mix them.
- **Write for the next person.** Assume they have not read the conversation you had, only the files.
- **Be specific.** Give numbers, settings and file names, not "it looked fine".
- **Record failures.** A dead end, written down, saves the next person a day.
- **Keep the cosmic language:** `mu`, `G`, `Re`, `Com`, `Sh`, `Ht`, `abs`, `E_cal`, `B_cal`.

---

## 7. Communication and rhythm

- **Asynchronous by default.** Issues for questions, pull requests for work, `team/LOG.md` for the journal.
- **Weekly summary.** Each Sunday the scribe writes a short entry at the top of `team/LOG.md`: what was done, what was confirmed, what is blocked, who needs help.
- **Disagreements** about results are settled by reproduction, not argument. Disagreements about direction are settled by the maintainer, recorded in `team/DECISIONS.md`.
- **Urgent problems** (CI red, a verified formula seems wrong, a rule is broken): open an issue with the label `urgent`.

---

## 8. Quick checklist

**Before starting:**

- [ ] read the reading list (section 4);
- [ ] verification scripts pass locally;
- [ ] scenario claimed (issue and board);
- [ ] notebook opened or continued.

**Before stopping a session:**

- [ ] notebook updated;
- [ ] log entry written;
- [ ] work pushed to your branch.

**Before delivering:**

- [ ] script, settings stamp, raw output and convergence check committed;
- [ ] T0 status stated;
- [ ] every statement labelled;
- [ ] findings index and board updated;
- [ ] "Scenario result" issue opened.
