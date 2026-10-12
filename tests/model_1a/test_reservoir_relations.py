import math

import pytest

from concentric_tube_cht.model_1a.reservoir_relations import (
    reservoir_to_inlet_state,
)


P0 = 10.0e5
T0 = 1000.0
GAMMA = 1.2105
GAS_CONSTANT = 268.10
DIAMETER = 0.010
AREA = math.pi * DIAMETER**2 / 4.0


def test_reservoir_to_inlet_reference_state():
    state = reservoir_to_inlet_state(
        M_in=0.5,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        area=AREA,
    )

    assert state.temperature == pytest.approx(
        974.3620973,
        rel=1.0e-9,
    )
    assert state.pressure == pytest.approx(
        861262.4162,
        rel=1.0e-9,
    )
    assert state.velocity == pytest.approx(
        281.1648324,
        rel=1.0e-9,
    )
    assert state.mass_flow == pytest.approx(
        0.07280633685,
        rel=1.0e-9,
    )


def test_reservoir_conditions_are_recovered():
    state = reservoir_to_inlet_state(
        M_in=0.5,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        area=AREA,
    )

    stagnation_factor = (
        1.0
        + 0.5 * (GAMMA - 1.0) * state.mach**2
    )

    reconstructed_T0 = state.temperature * stagnation_factor

    reconstructed_p0 = (
        state.pressure
        * stagnation_factor ** (GAMMA / (GAMMA - 1.0))
    )

    assert reconstructed_T0 == pytest.approx(T0, rel=1.0e-12)
    assert reconstructed_p0 == pytest.approx(P0, rel=1.0e-12)


def test_total_enthalpy_is_conserved():
    state = reservoir_to_inlet_state(
        M_in=0.5,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        area=AREA,
    )

    cp = GAMMA * GAS_CONSTANT / (GAMMA - 1.0)

    inlet_total_enthalpy = (
        cp * state.temperature
        + 0.5 * state.velocity**2
    )

    reservoir_enthalpy = cp * T0

    assert inlet_total_enthalpy == pytest.approx(
        reservoir_enthalpy,
        rel=1.0e-12,
    )