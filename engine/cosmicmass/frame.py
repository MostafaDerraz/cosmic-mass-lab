"""The leaning frame of the General Number plane (task E1.1).

A frame is the Composer Unit `G`: its length `s = |G|` and its frame angle
`phi_bar`. Compendium, section "General Numbers", paragraph "Rule of the Square
and the product": `G^2 = 2 s cos(phi_bar) G - s^2`.

The angle is stored in radians and read and written in degrees. `s` and
`phi_bar` are floats in normal use; in exact mode they are `mpmath.mpf` values,
and every operation then runs at the working precision of mpmath.
"""

import math
from dataclasses import dataclass

import mpmath

from . import constants as _constants

__all__ = ["Frame", "DEFAULT_PHI_BAR_DEG"]

# Frame angle held by the theory (compendium, "The cosmic mass and its unit").
DEFAULT_PHI_BAR_DEG = 79.86


def _is_exact(x) -> bool:
    return isinstance(x, mpmath.mpf)


@dataclass(frozen=True, slots=True)
class Frame:
    """`Frame(s, phi_bar, units)`: the Composer Unit, its length and its frame angle."""

    s: float  # |G| = sqrt(K/G_N) in SI, 1 in scaled units
    phi_bar: float  # frame angle in radians, strictly between 0 and pi
    units: str = "SI"  # "SI" or "scaled"

    def __post_init__(self):
        if self.units not in ("SI", "scaled"):
            raise ValueError(f"units must be 'SI' or 'scaled', not {self.units!r}")
        if not self.s > 0:
            raise ValueError(f"the length s must be positive, not {self.s!r}")
        if not 0 < self.phi_bar < mpmath.pi:
            raise ValueError("the frame angle must lie strictly between 0 and 180 degrees")

    # --- constructors -------------------------------------------------
    @classmethod
    def from_degrees(cls, s, phi_bar_deg, units: str = "SI") -> "Frame":
        """Build a frame from an angle in degrees. An `mpf` length gives an exact frame."""
        if _is_exact(s) or _is_exact(phi_bar_deg):
            return cls(mpmath.mpf(s), mpmath.radians(mpmath.mpf(phi_bar_deg)), units)
        return cls(float(s), math.radians(phi_bar_deg), units)

    @classmethod
    def si(cls, phi_bar_deg: float = DEFAULT_PHI_BAR_DEG) -> "Frame":
        """The SI frame: `s = sqrt(K/G_N) = 1.16043e10 kg/C`."""
        return cls.from_degrees(_constants.SI.s, phi_bar_deg, "SI")

    @classmethod
    def scaled(cls, phi_bar_deg: float = 90.0) -> "Frame":
        """The scaled frame: `s = 1`; the natural frame angle is 90 degrees."""
        return cls.from_degrees(_constants.SCALED.s, phi_bar_deg, "scaled")

    # --- the angle ----------------------------------------------------
    @property
    def exact(self) -> bool:
        """True when the frame holds mpmath values (exact mode)."""
        return _is_exact(self.s) or _is_exact(self.phi_bar)

    @property
    def phi_bar_deg(self):
        """The frame angle in degrees."""
        return mpmath.degrees(self.phi_bar) if self.exact else math.degrees(self.phi_bar)

    def cos_phi(self):
        return mpmath.cos(self.phi_bar) if self.exact else math.cos(self.phi_bar)

    def sin_phi(self):
        return mpmath.sin(self.phi_bar) if self.exact else math.sin(self.phi_bar)

    def square_rule(self) -> tuple:
        """Rule of the Square as a pair `(p, q)` with `G^2 = p + q G`:
        `p = -s^2` and `q = 2 s cos(phi_bar)`."""
        return (-self.s * self.s, 2 * self.s * self.cos_phi())

    # --- reading and writing (degrees) ---------------------------------
    def to_dict(self) -> dict:
        """Settings form of a float frame; the angle is written in degrees."""
        return {"s": float(self.s), "phi_bar_deg": float(self.phi_bar_deg), "units": self.units}

    @classmethod
    def from_dict(cls, data: dict) -> "Frame":
        return cls.from_degrees(data["s"], data["phi_bar_deg"], data.get("units", "SI"))
