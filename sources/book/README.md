# The Cosmic Mass, Volume 1: One Number for Mass and Charge (first draft)

LaTeX sources of the book by Mostafa Derraz, typeset with the IMISE classicthesis template
(https://github.com/IMISE/imise-classicthesis), which is based on "A Classic Thesis Style"
by André Miede and Ivo Pletikosić (GNU GPL, version 2 or later; see the header of classicthesis.sty).

## Build

    latexmk -pdf book.tex

pdfLaTeX and BibTeX, with a standard TeX Live or MiKTeX installation (KOMA-Script, tcolorbox, pgfplots, natbib, cleveref, longtable).

The bibliography is made by BibTeX from `book.bib`. If your editor runs only pdfLaTeX, run BibTeX once (or use `latexmk`), then pdfLaTeX twice; otherwise the citations print as question marks. A ready `book.bbl` is included, so a plain pdfLaTeX run already prints the bibliography until a new citation is added.

## Files

- `book.tex`: main file (three parts, chapters, appendices; Part IV, with Chapter 10, was added on 8 November 2026, and Chapter 11 on 9 November 2026)
- `classicthesis-config.tex`: layout options, book data (title, volume, author), packages, theorem styles, notation macros, the "In plain words" and "Worked example" boxes
- `classicthesis.sty`: the template style
- `book.bib`: bibliography
- `FrontBackmatter/`: title page, copyright page, epigraph, preface, how to read, contents, bibliography, colophon
- `Chapters/`: the chapters, the Coda and Appendices A to F. `Chapters/Chapter03-Algebra.tex` holds Section 4.9, which `Chapter03.tex` reads in at the end of Chapter 4. `Chapters/Chapter06/`, `Chapters/Chapter07/` and `Chapters/Chapter08/` hold the figures and tables of the chapters in those files, and `AppendixE-table.tex` the table of Appendix E
- `book-updated.pdf`: the whole book after the update of 6 November 2026 (285 pages); replaced on 8 November 2026 by the build with Chapter 10 (309 pages), and on 9 November 2026 by the build with Chapter 11 (325 pages), then by the build with Chapters 10 and 11 shortened (311 pages)
- `removed-in-shortening-2026-11-09/`: the full Chapters 10 and 11 as they stood before they were shortened on 9 November 2026, kept for reference and not compiled, with `CHANGES.md`, the list of what was removed
- `removed-in-chain-cleanup/`: the text taken out of Chapters 3 and 6 to 9 in the chain cleanup, one file per chapter, kept for reference and not compiled, with `CHANGES.md`, the list of every change
- `backup-before-2026-11-02-update/`, `backup-before-2026-11-03-update/`, `backup-before-2026-11-06-update/`: the files as they were before each update (the last one before the chain cleanup and Section 4.9)

## Structure, and which file holds which chapter

The file names keep the numbering of the first draft, so that nothing had to be renamed or deleted when the new Chapter 2 was added.

| Chapter | Title | File |
|---|---|---|
| Part I | Beyond the Right Angle | |
| 1 | One Number for Two Laws | `Chapters/Chapter01.tex` |
| 2 | Old Laws, New Results | `Chapters/Alignment.tex` |
| 3 | The Oblique Triangle | `Chapters/Chapter02.tex` |
| 4 | General Numbers | `Chapters/Chapter03.tex` and `Chapters/Chapter03-Algebra.tex` (Section 4.9) |
| Part II | The Foundations of Cosmic Mass Theory | |
| 5 | Prologue: Limits of Physics Language | `Chapters/Chapter04.tex` |
| 6 | The Interaction Law | `Chapters/Chapter05.tex` |
| 7 | The Cosmic Mass Standard Model | `Chapters/Chapter06.tex` and `Chapters/Chapter06/` |
| Part III | From Nuclei to Fields | |
| 8 | The Closure Law | `Chapters/Chapter07.tex` and `Chapters/Chapter07/` |
| 9 | The Unified Cosmic Law | `Chapters/Chapter08.tex` and `Chapters/Chapter08/` |
| Part IV | Open Questions | |
| 10 | What Is the Coupling Force? An Open Discussion | `Chapters/Chapter09.tex` |
| 11 | What the Frame Allows: Antimatter, Negative Mass and Black Holes | `Chapters/Chapter10.tex` |
| | Coda: The Next Rung | `Chapters/Coda.tex` |
| | Appendices A to F | `Chapters/AppendixA.tex` to `AppendixF.tex` |

This is Volume 1. It shows that the cosmic mass language can speak the old laws, all of them in one form; Volume 2 is for new predictions, beginning with the Coupling Force (Chapter 2, Section 2.4). Chapter 10 gathers the open readings of the Coupling Force, with their status and the questions they hand to Volume 2. Chapter 11 draws the whole frame (antimatter, the half of negative mass and the black holes), with what has been measured and what is only allowed.

## Notes for the author

- Placeholders are printed in red in square brackets. Search the sources for `\placeholder{`.
- Every formula of the manuscripts is reproduced exactly; Chapter 8 refits the error law of the Closure Law on measured masses, and says so.
- Chapters 4 and 6 keep the order and numbering of the manuscript, with the chapter number raised by one since the new Chapter 2: Axiom 3.1, Theorem 5.6 and Remark 5.15 of the papers are Axiom 4.1, Theorem 6.6 and Remark 6.15 here. The statements added after them continue the numbering (Definition 4.11 to Remark 4.60 in Chapter 4, Proposition 6.16 and Definition 6.17 in Chapter 6); Proposition 6.18 was added after them on 6 November 2026.
- Chapters 6 to 9 read every measured force, energy and momentum by its real part (Definition 6.17) and hold the angle of the Composer axis at 79.86 degrees, the angle at which the square of the Composer Unit points backward and like charges push apart (Sections 6.3.2 and 6.6.3). They discuss no other angle; Chapter 4 keeps the angle free as mathematics, and Section 7.2.4 shows how the angle shapes the ellipse.
- Chain cleanup (3 November 2026): Chapters 1, 2 and 5 are unchanged. In Chapters 3, 4 and 6 to 9 every part feeds a later chapter. Chapter 4 now holds all the mathematics that Chapters 6 to 9 use, including the Rule of the Square and the real part of a product (Section 4.6.3), and Table 4.1 lists each tool and where it is used.
- The algebra of General Numbers (6 November 2026): Section 4.9, the last section of Chapter 4, gives each law of the complex numbers in Appendix C its form in a relative plane, with its complex case at the angle pi/2, and Subsection 4.9.7 reads each result by its two parts and by its shadow and height. Appendix F lists the same laws in tables beside the complex ones. The physical chapters do not use this section; it stands beside the road to Part II. Its new commands (`\Jp`, `\Ht`, `\incomplex`, `\gnkeep`, and the table environment of Appendix F) are defined in the two new files with `\providecommand` and `\newenvironment`, so `classicthesis-config.tex` did not change for it.
- Readings from Section 4.9 (6 November 2026): the force read from two identity cards (Proposition 6.18, Section 6.5.4); what the Composer coordinate of a product holds, the Coupling term beside a much larger share of the term of two charges (Section 6.5.2, Section 7.6.2, Chapter 9 Step 1, the Coda); one cause for both failures of the shadow (Section 6.5.3); the frame as a chart of pulls and pushes (Section 7.6, reading 8); the Coulomb barrier read from the cards (Chapter 8, Step 2); bodies on one ray move alike (Chapter 9, Step 4). No formula of the manuscripts changed.
- Chapter 4 defines the words of the General Numbers in mathematics, and Chapters 6 to 9 use the same words in physics; Table 6.1 pairs them.
- Chapters 8 and 9 use the statement type `Law`, numbered with the theorems; the three lines at the top of each of the two files define it.
- `Chapters/Chapter06/fig-sign.tex` and `fig-sign-real.tex`, the two sign figures of earlier versions of Chapter 7, are kept but no longer used. `Chapters/Chapter03_old.tex` is not included in the book.
- Chapter 10 (8 November 2026): an open discussion of the Coupling Force, in the new Part IV. It states the author's two readings as Views, a statement type defined in the first lines of `Chapters/Chapter09.tex` (the Coupling term as two masses, and a particle as a product of particles), and adds the diagonal of the frame and its place in general relativity, the charge share, the flip of every charge, a body in a field, a table of the measurements that already limit each reading, and the questions for Volume 2. Each reading has a status: what holds as mathematics, what measurement excludes, what stays open. The chapter holds the angle at 79.86 degrees; it notes that the factors of the General Sinus Transform have angle sums (79.86 + 180 k)/N degrees, and it asks, as an open question and without giving any other angle, whether 79.86 degrees belongs to nature or to the units. Two labels were added to the first two products of Section A.9, so that Chapter 10 can cite equations (A.64) and (A.65); no formula changed. The Coda and the Preface gained one sentence each, How to Read two (one for Part IV, one for the Views), and `book.bib` seventeen references.
- Chapter 11 (9 November 2026): what the frame allows, in seven steps, in Part IV after Chapter 10. The antiparticle of a body has the cosmic mass m - qG, the skew flip of Remark 4.39 and not the conjugate, so on the frame it is the mirror image of its particle in the mass axis. The chapter reads a particle and its antiparticle together (the pull, annihilation, pair creation); shows that two numbers cannot tell neutral matter from neutral antimatter; reads by the law of motion of Chapter 9 what a body of negative mass would do (the runaway pair that Bondi found in Einstein's theory); and writes the horizons of a charged black hole, a result of Einstein's theory that it quotes and does not derive, with the real part of the square of the cosmic mass. Every place of the frame gets a status (measured, allowed, or empty so far), and four questions go to Volume 2. The chapter holds the angle at 79.86 degrees and discusses no other angle. The Preface and the Coda gained one sentence each, the Part IV sentence of How to Read now names Chapter 11 as well, the last sentence of Chapter 10 points to Chapter 11, and `book.bib` gained four references. No earlier formula changed.
- Shortening (9 November 2026): at the author's request, Chapters 10 and 11 were made shorter and more direct, Chapter 10 from 20 pages to 12 and Chapter 11 from 16 to 9. Every View, definition, proposition, remark and proof is kept word for word, under the same numbers; the proof of Proposition 10.7 now quotes the classical product of sines instead of deriving it. Removed are the road maps, the "In plain words" boxes inside the steps, the "Looking Back" sections, Figures 10.2 and 11.2, Table 11.1 and the longer worked examples; each chapter keeps one closing box. The full earlier chapters are in `removed-in-shortening-2026-11-09/`, with `CHANGES.md`.
- Citing the General Numbers paper (9 November 2026): the paper on General Numbers is published on Zenodo (DOI 10.5281/zenodo.24091305; key `DerrazGNPreprint2026` in `book.bib`). The book cites it in seven places: the Sources on the copyright page, the Preface, Section 4.1, the pi/5 worked example, the spiral of powers, Section 4.9 and Appendix F. Nothing was removed and no formula changed.
- Sources cleaned (9 November 2026): at the author's request, the book cites only two of the author's works, General Sinus Principles (Research Square) and the General Numbers paper (Zenodo). The other works (the General Numbers manuscript, the Royal Society manuscript, Stations and Meditations, the Closure Law draft, the fusion guide, the Unified Cosmic Law draft and the 2017 General Sinus article) are drafts of this book: they are neither cited nor named. Their entries stay in `book.bib`, unused. No formula changed, and the workbooks stay.
- `drafting=true` in `classicthesis-config.tex` prints the version line at the foot of each page and marks overfull lines. Set it to `false` for the final book.
