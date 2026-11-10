"""Physical constants of the virtual universe (task E1.1).

Source of the SI values: `verification/cosmic-mass/constants.py` (CODATA 2022),
the set against which the Cosmic Mass paper is verified. Formulas: compendium,
section "Cosmic Mass", paragraphs "The cosmic mass and its unit" (`s = sqrt(K/G)`)
and "Stoney mass" (`m_S = e s`).

Two sets exist (decision D7): SI, the reference, and scaled units, where
`G_N = K = 1` and therefore `s = 1`.
"""

import math
from dataclasses import dataclass

__all__ = ["Constants", "SI", "SCALED", "scaled"]

_EPS0 = 8.8541878188e-12  # vacuum electric permittivity, F/m (CODATA 2022)


@dataclass(frozen=True, slots=True)
class Constants:
    """One consistent set of constants. `s` and the Stoney mass are derived, never stored."""

    G_N: float  # Newtonian constant of gravitation
    K: float  # Coulomb constant
    c: float  # speed of light
    e: float  # elementary charge
    units: str  # "SI" or "scaled"

    def __post_init__(self):
        if self.units not in ("SI", "scaled"):
            raise ValueError(f"units must be 'SI' or 'scaled', not {self.units!r}")
        if min(self.G_N, self.K, self.c, self.e) <= 0:
            raise ValueError("G_N, K, c and e must be positive")

    @property
    def s(self) -> float:
        """Length of the Composer Unit, `s = |G| = sqrt(K/G_N)` (kg/C in SI)."""
        return math.sqrt(self.K / self.G_N)

    @property
    def stoney_mass(self) -> float:
        """Stoney mass `m_S = e s` (kg in SI)."""
        return self.e * self.s


SI = Constants(
    G_N=6.67430e-11,  # m^3 kg^-1 s^-2
    K=1.0 / (4.0 * math.pi * _EPS0),  # N m^2 C^-2
    c=299792458.0,  # m/s (exact)
    e=1.602176634e-19,  # C (exact)
    units="SI",
)


def scaled(c: float = 1.0, e: float = 1.0) -> Constants:
    """Scaled units: `G_N = K = 1`, so `s = 1` (D7).

    D7 fixes only `G_N` and `K`. The values of `c` and `e` in scaled units are a
    choice of the scenario; the default `c = e = 1` (Stoney units) is a labelled
    assumption, recorded in `team/notes/E1.md`.
    """
    return Constants(G_N=1.0, K=1.0, c=c, e=e, units="scaled")


SCALED = scaled()
