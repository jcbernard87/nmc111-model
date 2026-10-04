"""The protocol driver's constant-voltage step, with stand-in steppers (no model needed)."""
import pytest

from nmc_model.driver import SolverFailure, cv_step


class _P:
    i_1C = 1.0


class Jump:
    """V(I) jumps over the set voltage at I = 0.5: no current holds V = 0.45."""
    p = _P()

    def newton_step(self, state, h, I):
        return I

    def voltage(self, state, I):
        return 1.0 - I if I < 0.5 else 0.4 - I


def test_cv_never_accepts_a_state_off_the_set_voltage():
    """When the bracket collapses on the jump, the last state is 0.05 V off: a failure, not a CV state (#1)."""
    with pytest.raises(SolverFailure):
        cv_step(Jump(), None, 10.0, 0.45, 0.0)
