"""The readings of a General Number (task E1.4).

For `z = a + b G` in a frame `(s, phi_bar)` (compendium, section "General
Numbers", paragraph "Readings"; `ARCHITECTURE.md`, section 5):

- `Re(z) = a`, the Real coordinate;
- `Com(z) = b`, the Composer coordinate (the compendium writes the Composer part
  as `b G`; the engine returns its coordinate `b`);
- `Sh(z) = a + b s cos(phi_bar)`, the shadow;
- `Ht(z) = b s sin(phi_bar)`, the height;
- `abs(z) = sqrt(Sh^2 + Ht^2)`, the length, with `|z|^2 = a^2 + 2ab s cos(phi_bar) + b^2 s^2`.

Every function accepts a `GeneralNumber` (floats or exact mpmath values) or a
`GNArray`, and returns a real number or a real array.
"""

import math

import mpmath
import numpy as np

__all__ = ["Re", "Com", "Sh", "Ht", "abs"]


def Re(z):
    """The Real coordinate `a`."""
    return z.re


def Com(z):
    """The Composer coordinate `b`."""
    return z.com


def Sh(z):
    """The shadow `a + b s cos(phi_bar)`."""
    return z.re + z.com * (z.frame.s * z.frame.cos_phi())


def Ht(z):
    """The height `b s sin(phi_bar)`."""
    return z.com * (z.frame.s * z.frame.sin_phi())


def abs(z):  # noqa: A001  (the name of the reading in the theory)
    """The length `sqrt(Sh^2 + Ht^2)`."""
    sh, ht = Sh(z), Ht(z)
    if isinstance(sh, np.ndarray):
        return np.hypot(sh, ht)
    if z.frame.exact:
        return mpmath.sqrt(sh * sh + ht * ht)
    return math.hypot(sh, ht)
