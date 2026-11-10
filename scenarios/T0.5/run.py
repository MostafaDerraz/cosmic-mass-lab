"""T0.5: Field split: the Composer part of one proton's field is Coulomb's field.

Setup and pass criterion: docs/05-scenarios-and-tests.md, section 2, T0. Pre-registration:
team/notes/T0.md. The setup is built in engine/cosmicmass/gates/t0.py and stamped per case.

    python scenarios/T0.5/run.py        # writes results/T0.5/<run-id>/
"""

import sys

from cosmicmass.gates.t0 import t0_5, write_gate_result

PRODUCER = "Claude chat session on claude.ai / Claude Fable 5.2 (sponsor @MostafaDerraz)"

if __name__ == "__main__":
    outcome = t0_5()
    folder = write_gate_result(outcome, producer=sys.argv[1] if len(sys.argv) > 1 else PRODUCER)
    for name, case in outcome["cases"].items():
        print(f"{'pass' if case['pass'] else 'FAIL'}  T0.5 {name}")
    print(f"T0.5: {'pass' if outcome['pass'] else 'FAIL'}; results in {folder}")
    sys.exit(0 if outcome["pass"] else 1)
