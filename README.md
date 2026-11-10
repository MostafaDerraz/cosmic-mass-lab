# Cosmic Mass Lab

Working repository of Mostafa Derraz for the General Sinus, General Numbers and Cosmic Mass theory, and for the virtual-universe simulator built on them.

| Folder | Content |
|---|---|
| `sources/book/` | *The Cosmic Mass, Volume 1*: LaTeX sources (build with `latexmk -pdf book.tex`) and the built PDF |
| `sources/cosmic-mass-paper/` | Cosmic Mass paper, Zenodo version, doi:10.5281/zenodo.24094855 |
| `sources/general-numbers-paper/` | General Numbers paper, Zenodo version, doi:10.5281/zenodo.24091305, and the Gazette version |
| `sources/general-sinus-paper/` | General Sinus, revision 3 |
| `verification/` | verification scripts for every source, with their latest outputs |
| `docs/` | report on verified formulas, formulas compendium (LaTeX and PDF), scenarios and tests, engine architecture (Mermaid), general note |
| `simulator/` | the virtual universe of cosmic masses (to be built; see `docs/05-scenarios-and-tests.md`) |

## Verification

Every push runs all verification scripts (GitHub Actions, `.github/workflows/verify.yml`).

| Script | Checks |
|---|---|
| `verification/cosmic-mass/verify_cmt.py` | 106 |
| `verification/general-numbers/verify_paper.py` | 204 |
| `verification/general-numbers/book-algebra/verify_general_number_algebra.py` | 253 |
| `verification/general-numbers/book-laws-test/test_general_number_laws.py` | 144 laws |
| `verification/general-numbers/gazette-version/check_gazette.py` | 9,304 |
| `verification/general-sinus/verify_gsp.py` | 79 |

Run them locally with Python 3 and SymPy (`pip install sympy mpmath`), using the commands in the workflow file.

## Status

This repository is private. The papers are under review or posted as preprints; do not redistribute without the author's permission.
