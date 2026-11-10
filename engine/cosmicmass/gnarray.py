"""General Numbers on NumPy arrays (task E1.3).

`GNArray(re, com, frame)` holds many General Numbers `re + com G` of one frame as
two float arrays of the same shape. It is the type of the hot path: scalar fields
have shape `(N,)` and vector fields `(N, d)`. The formulas and the order of the
floating-point operations are those of `gn.py` (compendium, section "General
Numbers"), so an array operation gives, element by element, the same bits as the
scalar class.

Shapes follow NumPy broadcasting. To multiply N numbers by an `(N, d)` vector
field, index the numbers with `[:, None]`.
"""

import numpy as np

from .frame import Frame
from .gn import GeneralNumber

__all__ = ["GNArray"]


class GNArray:
    """Arrays of General Numbers `re + com G` in the frame `frame` (floats only)."""

    __slots__ = ("re", "com", "frame")
    __array_ufunc__ = None  # `ndarray * GNArray` is handled here, not by NumPy
    __hash__ = None

    def __init__(self, re, com, frame: Frame):
        if frame.exact:
            raise TypeError("GNArray needs a float frame; exact mode is for GeneralNumber")
        re, com = np.asarray(re, dtype=float), np.asarray(com, dtype=float)
        if re.shape != com.shape:
            re, com = (np.array(x) for x in np.broadcast_arrays(re, com))
        self.re, self.com, self.frame = re, com, frame

    # --- constructors and access ---------------------------------------
    @classmethod
    def zeros(cls, shape, frame: Frame) -> "GNArray":
        return cls(np.zeros(shape), np.zeros(shape), frame)

    @classmethod
    def from_scalars(cls, numbers) -> "GNArray":
        """A one-dimensional array from a non-empty sequence of `GeneralNumber` of one frame."""
        numbers = list(numbers)
        out = cls([z.re for z in numbers], [z.com for z in numbers], numbers[0].frame)
        if any(z.frame != out.frame for z in numbers):
            raise ValueError("General Numbers of different frames cannot be combined")
        return out

    @property
    def shape(self) -> tuple:
        return self.re.shape

    def __len__(self) -> int:
        return len(self.re)

    def __getitem__(self, index) -> "GNArray":
        return GNArray(self.re[index], self.com[index], self.frame)

    def item(self, *index) -> GeneralNumber:
        """One element as a scalar `GeneralNumber`."""
        return GeneralNumber(float(self.re[index]), float(self.com[index]), self.frame)

    def _coerce(self, other):
        """`(re, com)` of `other` in this frame, or None for an unknown type. Real
        scalars and real arrays are General Numbers with a zero Composer coordinate."""
        if isinstance(other, (GNArray, GeneralNumber)):
            if other.frame != self.frame:
                raise ValueError("General Numbers of different frames cannot be combined")
            return other.re, other.com
        if isinstance(other, (bool, np.bool_)):
            return None
        if isinstance(other, (int, float, np.integer, np.floating)) or (
            isinstance(other, np.ndarray) and other.dtype.kind in "iuf"
        ):
            return other, 0.0
        return None

    # --- sum, difference ---------------------------------------------------
    def add(self, other) -> "GNArray":
        w = self._coerce(other)
        if w is None:
            return NotImplemented
        return GNArray(self.re + w[0], self.com + w[1], self.frame)

    def __neg__(self) -> "GNArray":
        return GNArray(-self.re, -self.com, self.frame)

    def __sub__(self, other):
        w = self._coerce(other)
        if w is None:
            return NotImplemented
        return GNArray(self.re - w[0], self.com - w[1], self.frame)

    def __rsub__(self, other):
        return (-self).add(other)

    # --- product ----------------------------------------------------------
    def mul(self, other) -> "GNArray":
        """Product with the Rule of the Square, element by element (both coordinates)."""
        w = self._coerce(other)
        if w is None:
            return NotImplemented
        a, b, (c, d) = self.re, self.com, w
        p, q = self.frame.square_rule()  # G^2 = p + q G
        bd = b * d
        return GNArray(a * c + p * bd, (a * d + b * c) + q * bd, self.frame)

    # --- conjugate, length, inverse, division --------------------------------
    def conjugate(self) -> "GNArray":
        _, q = self.frame.square_rule()
        return GNArray(self.re + q * self.com, -self.com, self.frame)

    def norm_sq(self) -> np.ndarray:
        """`|z|^2 = a^2 + 2ab s cos(phi_bar) + b^2 s^2`, element by element."""
        p, q = self.frame.square_rule()
        return self.re * self.re + q * self.re * self.com - p * self.com * self.com

    def inverse(self) -> "GNArray":
        n = self.norm_sq()
        if np.any(n == 0):
            raise ZeroDivisionError("the General Number zero has no inverse")
        zbar = self.conjugate()
        return GNArray(zbar.re / n, zbar.com / n, self.frame)

    def __truediv__(self, other):
        w = self._coerce(other)
        if w is None:
            return NotImplemented
        return self.mul(GNArray(w[0], w[1], self.frame).inverse())

    def __rtruediv__(self, other):
        return self.inverse().mul(other)

    __add__ = __radd__ = add
    __mul__ = __rmul__ = mul

    # --- comparison ---------------------------------------------------------
    def __eq__(self, other):
        """Element-by-element exact equality (a Boolean array)."""
        try:
            w = self._coerce(other)
        except ValueError:
            return False
        if w is None:
            return NotImplemented
        return (self.re == w[0]) & (self.com == w[1])

    def allclose(self, other, rel_tol: float = 1e-12, abs_tol: float = 0.0) -> bool:
        """True when every element is close: the length of the difference is at most
        `max(rel_tol * max(|z|, |w|), abs_tol)`, as in `GeneralNumber.isclose`."""
        w = self._coerce(other)
        if w is None:
            raise TypeError(f"cannot compare a GNArray with {type(other).__name__}")
        w = GNArray(w[0], w[1], self.frame)
        size = lambda z: np.sqrt(np.maximum(z.norm_sq(), 0.0))  # noqa: E731
        bound = np.maximum(rel_tol * np.maximum(size(self), size(w)), abs_tol)
        return bool(np.all(size(self - w) <= bound))

    def __repr__(self):
        return f"GNArray(shape={self.shape}, phi_bar={self.frame.phi_bar_deg:g} deg, s={self.frame.s:g})"
