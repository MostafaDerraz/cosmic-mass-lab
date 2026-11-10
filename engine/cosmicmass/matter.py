"""Cosmic masses and the matter presets (task E1.5).

A body carries one number, its cosmic mass `mu = m + q G` (compendium, section
"Cosmic Mass", paragraph "The cosmic mass and its unit"). The presets are the
table "Matter presets" of `docs/05-scenarios-and-tests.md`, section 0, in SI
units (kg and C), with the values of that table.
"""

from dataclasses import dataclass, replace

from .frame import Frame
from .gn import GeneralNumber

__all__ = ["CosmicMass", "PRESETS", "preset"]


@dataclass(frozen=True, slots=True)
class CosmicMass:
    """One body: mass `m` (kg), charge `q` (C) and an optional name."""

    m: float
    q: float
    name: str = ""

    def mu(self, frame: Frame) -> GeneralNumber:
        """The cosmic mass `mu = m + q G` as a General Number of `frame`."""
        return GeneralNumber(self.m, self.q, frame)

    def with_charge(self, q: float) -> "CosmicMass":
        """The same body with another charge (for example a charged planet)."""
        return replace(self, q=q)


_E = 1.602176634e-19  # elementary charge, C

PRESETS = {
    body.name: body
    for body in (
        CosmicMass(9.1093837e-31, -_E, "electron"),
        CosmicMass(1.67262192e-27, +_E, "proton"),
        CosmicMass(1.67492750e-27, 0.0, "neutron"),
        CosmicMass(6.6446573e-27, 3.204353268e-19, "alpha"),
        CosmicMass(1e-15, +1e-16, "dust grain +"),
        CosmicMass(1e-15, -1e-16, "dust grain -"),
        CosmicMass(7.342e22, 0.0, "moon"),
        CosmicMass(5.972e24, 0.0, "Earth"),
        CosmicMass(1.989e30, 0.0, "Sun"),
        CosmicMass(5.683e26, 0.0, "Saturn"),
    )
}


def preset(name: str) -> CosmicMass:
    """A matter preset by name; the names are the keys of `PRESETS`."""
    try:
        return PRESETS[name]
    except KeyError:
        raise KeyError(f"unknown matter preset {name!r}; known: {sorted(PRESETS)}") from None
