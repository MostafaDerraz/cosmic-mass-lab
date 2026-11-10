"""General Numbers: scalar arithmetic in a leaning frame (task E1.2).

A General Number is `z = a + b G`, with `a` its Real coordinate and `b` its
Composer coordinate, in a `Frame(s, phi_bar)`. Every formula here is from the
compendium, section "General Numbers":

- product ("Rule of the Square and the product"):
  `(a + bG)(c + dG) = (ac - s^2 bd) + (ad + bc + 2 s cos(phi_bar) bd) G`
- squared length ("Readings"): `|z|^2 = a^2 + 2ab s cos(phi_bar) + b^2 s^2`
- conjugate and inverse ("Conjugate, skew flip, division"):
  `conj(z) = (a + 2b s cos(phi_bar)) - bG` and `1/z = conj(z) / |z|^2`

This scalar class is for tests, exact checks and readability. The hot path uses
arrays (`gnarray.py`, task E1.3). Exact mode: give `mpmath.mpf` coordinates and
an exact frame, and every operation runs at the working precision of mpmath.
"""

import math
from dataclasses import dataclass
from numbers import Real

import mpmath

from .frame import Frame

__all__ = ["GeneralNumber"]

_SCALARS = (Real, mpmath.mpf)


@dataclass(frozen=True, slots=True, eq=False)
class GeneralNumber:
    """`z = re + com G` in the frame `frame`."""

    re: float
    com: float
    frame: Frame

    # --- constructors -------------------------------------------------
    @classmethod
    def unit(cls, frame: Frame) -> "GeneralNumber":
        """The Composer Unit `G` itself: `0 + 1 G`."""
        return cls(0 * frame.s, frame.s / frame.s, frame)

    def _coerce(self, other):
        """Return `other` as a General Number of this frame, or None if it is neither
        a real scalar nor a General Number."""
        if isinstance(other, GeneralNumber):
            if other.frame != self.frame:
                raise ValueError("General Numbers of different frames cannot be combined")
            return other
        if isinstance(other, _SCALARS) and not isinstance(other, bool):
            return GeneralNumber(other, 0 * other, self.frame)
        return None

    # --- addition and subtraction ---------------------------------------
    def add(self, other) -> "GeneralNumber":
        w = self._coerce(other)
        if w is None:
            return NotImplemented
        return GeneralNumber(self.re + w.re, self.com + w.com, self.frame)

    def __neg__(self) -> "GeneralNumber":
        return GeneralNumber(-self.re, -self.com, self.frame)

    def __sub__(self, other):
        w = self._coerce(other)
        return NotImplemented if w is None else self.add(-w)

    def __rsub__(self, other):
        return (-self).add(other)

    # --- product --------------------------------------------------------
    def mul(self, other) -> "GeneralNumber":
        """Product with the Rule of the Square (both coordinates)."""
        w = self._coerce(other)
        if w is None:
            return NotImplemented
        a, b, c, d = self.re, self.com, w.re, w.com
        p, q = self.frame.square_rule()  # G^2 = p + q G
        bd = b * d  # computed once, so that z * w == w * z bit for bit
        return GeneralNumber(a * c + p * bd, (a * d + b * c) + q * bd, self.frame)

    # --- conjugate, length, inverse, division -------------------------------
    def conjugate(self) -> "GeneralNumber":
        """`conj(z) = (a + 2b s cos(phi_bar)) - bG`, so that `z conj(z) = |z|^2`."""
        _, q = self.frame.square_rule()
        return GeneralNumber(self.re + q * self.com, -self.com, self.frame)

    def norm_sq(self):
        """`|z|^2 = a^2 + 2ab s cos(phi_bar) + b^2 s^2` (a real number, zero only for z = 0)."""
        p, q = self.frame.square_rule()
        return self.re * self.re + q * self.re * self.com - p * self.com * self.com

    def inverse(self) -> "GeneralNumber":
        """`1/z = conj(z) / |z|^2`."""
        n = self.norm_sq()
        if n == 0:
            raise ZeroDivisionError("the General Number zero has no inverse")
        zbar = self.conjugate()
        return GeneralNumber(zbar.re / n, zbar.com / n, self.frame)

    def __truediv__(self, other):
        w = self._coerce(other)
        return NotImplemented if w is None else self.mul(w.inverse())

    def __rtruediv__(self, other):
        return self.inverse().mul(other)

    __add__ = __radd__ = add
    __mul__ = __rmul__ = mul

    # --- equality -----------------------------------------------------------
    def __eq__(self, other):
        """Exact equality of both coordinates (a real scalar `x` equals `x + 0G`)."""
        try:
            w = self._coerce(other)
        except ValueError:
            return False
        if w is None:
            return NotImplemented
        return self.re == w.re and self.com == w.com

    def __hash__(self):
        return hash((self.re, self.com, self.frame))

    def isclose(self, other, rel_tol: float = 1e-12, abs_tol: float = 0.0) -> bool:
        """Equality with tolerance: the length of the difference is at most
        `max(rel_tol * max(|z|, |w|), abs_tol)`. The length is the frame's own `abs`."""
        w = self._coerce(other)
        if w is None:
            raise TypeError(f"cannot compare a General Number with {type(other).__name__}")
        sqrt = mpmath.sqrt if self.frame.exact else math.sqrt
        # max(..., 0) guards the square root against a rounding error below zero
        size = lambda z: sqrt(max(z.norm_sq(), 0))  # noqa: E731
        return size(self - w) <= max(rel_tol * max(size(self), size(w)), abs_tol)

    def __repr__(self):
        return f"GeneralNumber({self.re!r}, {self.com!r}, phi_bar={float(self.frame.phi_bar_deg):g} deg, s={float(self.frame.s):g})"
