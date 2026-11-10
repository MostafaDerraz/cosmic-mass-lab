"""L0, the reference law (task E4.2; decision D1).

Compendium, section "Cosmic Mass", paragraph "Cosmic Field Law and law of motion":

    Re(mu) dv/dt = Re[mu (E_cal + v x B_cal)]
    <=>  m dv/dt = m (g + v x B_g) + q (E + v x B),   with E = -s^2 e and B = -s^2 b.

L0 reads only the Real coordinate of the product; the Composer coordinate is not
read as a force. The Real coordinate does not depend on the frame angle.
"""

import numpy as np

from ..fields.sample import FieldSample
from ..frame import Frame
from ..gnarray import GNArray
from ..state import State
from .base import ForcePair, Law, LawSpec, register

__all__ = ["L0", "lorentz_argument", "cross"]


def cross(v: np.ndarray, B: np.ndarray) -> np.ndarray:
    """`v x B` for N bodies. In 3D, `B` has shape `(N, 3)`. In 2D, `B` has shape `(N,)`
    and points out of the plane, so `v x B = (v_y B, -v_x B)`."""
    if v.shape[1] == 3:
        return np.cross(v, B)
    return np.stack([v[:, 1] * B, -v[:, 0] * B], axis=1)


def lorentz_argument(state: State, fields: FieldSample, frame: Frame) -> GNArray:
    """`E_cal + v x B_cal` at each body, a General Number vector field of shape `(N, d)`."""
    return GNArray(fields.g + cross(state.v, fields.B_g), fields.e + cross(state.v, fields.b), frame)


@register
class L0(Law):
    """The paper's law: force `Re[mu (E_cal + v x B_cal)]`, inertia `Re(mu) = m`."""

    reading = "L0"
    hypothesis = False

    def __init__(self, spec: LawSpec | None = None):
        super().__init__(spec)
        if not self.spec.is_default_except_reading():
            raise ValueError("L0 is the reference law: it takes no L2-L8 option")

    def forces(self, state: State, fields: FieldSample, frame: Frame) -> ForcePair:
        product = state.mu(frame)[:, None] * lorentz_argument(state, fields, frame)
        return ForcePair(real=product.re, composer=np.zeros_like(product.re))
