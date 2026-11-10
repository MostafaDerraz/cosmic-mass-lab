"""The environment of the bodies: walls (task E6a; reflecting walls only in the MVP).

A reflecting wall turns back the velocity component of a body that has crossed it
and is still moving outwards. The body is not moved, so its position, its speed and
therefore its kinetic and potential energy are exactly what they were: reflecting
walls conserve energy exactly. The momentum given to each wall is reported, for the
pressure on the walls (E8.1). Absorbing walls, wrap-around, inflow and the heat
bath (E6a, E6b) are outside the MVP.
"""

import numpy as np

from .state import State

__all__ = ["Environment"]


class Environment:
    """`boundary` is "open" or "reflecting". A reflecting box is given by its `lower`
    and `upper` corners, one number per dimension."""

    def __init__(self, boundary: str = "open", lower=None, upper=None):
        if boundary not in ("open", "reflecting"):
            raise ValueError("boundary must be 'open' or 'reflecting' in the MVP")
        self.boundary = boundary
        self.lower = self.upper = None
        if boundary == "reflecting":
            if lower is None or upper is None:
                raise ValueError("a reflecting box needs its 'lower' and 'upper' corners")
            self.lower, self.upper = np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)
            if self.lower.shape != self.upper.shape or self.lower.ndim != 1 or not np.all(self.lower < self.upper):
                raise ValueError("'lower' and 'upper' must be lists of one length, with lower < upper")

    def apply(self, state: State, dt: float, rng=None) -> list:
        """Apply the walls to `state`, in place. Returns one event per wall that was hit:
        `{"wall": (axis, side), "bodies": count, "impulse": momentum given to the wall}`."""
        events = []
        if self.boundary == "open":
            return events
        if self.lower.shape[0] != state.dimension:
            raise ValueError(f"the box has {self.lower.shape[0]} dimensions, the state has {state.dimension}")
        moving = state.alive & ~state.fixed
        for axis in range(state.dimension):
            for side, bound, sign in (("lower", self.lower[axis], -1.0), ("upper", self.upper[axis], 1.0)):
                out = moving & (sign * (state.x[:, axis] - bound) > 0) & (sign * state.v[:, axis] > 0)
                if np.any(out):
                    impulse = float(2.0 * np.sum(state.m[out] * state.v[out, axis]))
                    state.v[out, axis] = -state.v[out, axis]
                    events.append({"wall": (axis, side), "bodies": int(out.sum()), "impulse": impulse})
        return events
