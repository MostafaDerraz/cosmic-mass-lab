# General Sinus Principles: review and rebuild plan (10 November 2026)

## Errors found by numerical audit of rev1 (all other identities in Sections 4-5 pass)
1. `tan(2a,phi) - sec(2a,phi) = sin(a,phi)` (Remark after double-angle): FALSE. Delete.
2. `tan(a,phi) = sin(a/2,phi) + sec(a,phi)` (half-angle prop.): FALSE. Delete.
3. Second partials: `sin_yy`, `cos_xx` must be `f*(1 + 2cot^2(x+y))`, not `(2cot^2 - 1)`;
   `tan_yy = tan(x,y)*(1 + 2cot^2 y)` (paper uses x+y and -1). Fix.
4. Derivative along D = d/dx - d/dy (x+y fixed, the paper's "Axiom"):
   - paper's `sin^(2) = +sin`, `cos^(2) = +cos` and Remark "without sign change": FALSE; true is `-sin`, `-cos`.
   - paper's n-th derivative via omega(z,n) is wrong for n = 2, 3 (mod 4); n-th primitive likewise.
   - CORRECT, clean theorem: `D^n f(x,y) = f(x + n*pi/2, y - n*pi/2)` for f = sin, cos (all n, also negative n for primitives).
   - paper's `tan^(n) = n! tan' tan(pi/2+y,y)^(n-1)` holds only for n = 1, 2. Delete or replace.
5. Section 6 (product identities / GST) already fixed in rev1 and verified.

## Rebuild plan (rev2)
- Reframe the "Axiom" of opposite directions: D is the derivative with the oblique frame fixed
  (x+y = constant), i.e. moving the point around the circle. State it as a definition, not an axiom.
- Replace the omega operator by the shift theorem above (one line proof: D^n sin x = sin(x+n pi/2), sin(x+y) constant).
- Cut to ~16-18 pages: keep 3-4 figures (composed circle, oblique triangle, coordinate reading, addition configuration);
  drop value tables, 2-D/3-D surface galleries, the "common mistake" figure, the tangent/cotangent circles.
- Drop the example sentence "a value inaccessible to the relative function".
- Pick one sign in the tan addition / double formulas (both signs are true; the +- notation confuses).
- Add a verification script (identities, derivatives, products) and cite it.
- Target journal (honest): Elemente der Mathematik or Journal for Geometry and Graphics; or a short note
  (product identities) for The Mathematical Gazette after its decision on submission 262553062.
