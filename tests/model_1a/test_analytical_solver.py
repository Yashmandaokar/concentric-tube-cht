import math
import pytest

from concentric_tube_cht.model_1a.analytical_solver import (
    solve_critical_fanno_flow,
    solve_fanno_flow_for_back_pressure,
    solve_trial_fanno_flow,
)


P0 = 10.0e5
T0 = 1000.0
GAMMA = 1.2105
GAS_CONSTANT = 268.10
DIAMETER = 0.010
LENGTH = 0.020
FRICTION_FACTOR = 0.01430

AREA = math.pi * DIAMETER**2 / 4.0
CP = GAMMA * GAS_CONSTANT / (GAMMA - 1.0)


def make_trial_solution():
    return solve_trial_fanno_flow(
        M_in=0.5,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )


def test_trial_solution_remains_subsonic():
    solution = make_trial_solution()

    assert 0.5 < solution.exit_mach < 1.0
    assert not solution.choked_at_exit


def test_trial_solution_conserves_mass():
    solution = make_trial_solution()

    inlet_mass_flow = (
        solution.inlet.density
        * solution.inlet.velocity
        * AREA
    )

    exit_mass_flow = (
        solution.exit_density
        * solution.exit_velocity
        * AREA
    )

    assert inlet_mass_flow == pytest.approx(
        solution.inlet.mass_flow,
        rel=1.0e-12,
    )

    assert exit_mass_flow == pytest.approx(
        inlet_mass_flow,
        rel=1.0e-12,
    )


def test_trial_solution_conserves_total_enthalpy():
    solution = make_trial_solution()

    inlet_total_enthalpy = (
        CP * solution.inlet.temperature
        + 0.5 * solution.inlet.velocity**2
    )

    exit_total_enthalpy = (
        CP * solution.exit_temperature
        + 0.5 * solution.exit_velocity**2
    )

    assert exit_total_enthalpy == pytest.approx(
        inlet_total_enthalpy,
        rel=1.0e-12,
    )

    assert inlet_total_enthalpy == pytest.approx(
        CP * T0,
        rel=1.0e-12,
    )


def test_trial_solution_satisfies_perfect_gas_equation():
    solution = make_trial_solution()

    reconstructed_exit_pressure = (
        solution.exit_density
        * GAS_CONSTANT
        * solution.exit_temperature
    )

    assert reconstructed_exit_pressure == pytest.approx(
        solution.exit_pressure,
        rel=1.0e-12,
    )


def test_friction_reduces_static_and_stagnation_pressure():
    solution = make_trial_solution()

    assert solution.exit_pressure < solution.inlet.pressure
    assert solution.exit_stagnation_pressure < P0

def test_critical_solution_is_sonic_at_exit():
    solution = solve_critical_fanno_flow(
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert solution.inlet.mach == pytest.approx(
        0.876996273582931,
        rel=1.0e-10,
    )

    assert solution.exit_mach == pytest.approx(
        1.0,
        abs=1.0e-12,
    )

    assert solution.choked_at_exit


def test_critical_solution_reference_values():
    solution = solve_critical_fanno_flow(
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert solution.inlet.mass_flow == pytest.approx(
        0.09725721294,
        rel=1.0e-9,
    )

    assert solution.exit_temperature == pytest.approx(
        904.7726759,
        rel=1.0e-9,
    )

    assert solution.exit_pressure == pytest.approx(
        554329.1116,
        rel=1.0e-9,
    )

def test_atmospheric_back_pressure_selects_choked_branch():
    solution = solve_fanno_flow_for_back_pressure(
        p0=P0,
        T0=T0,
        back_pressure=1.013e5,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert solution.regime == "choked"
    assert solution.flow.choked_at_exit
    assert solution.flow.exit_mach == pytest.approx(1.0)
    assert solution.flow.exit_pressure > solution.back_pressure


def test_unchoked_branch_matches_back_pressure():
    solution = solve_fanno_flow_for_back_pressure(
        p0=P0,
        T0=T0,
        back_pressure=8.0e5,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert solution.regime == "unchoked"
    assert not solution.flow.choked_at_exit
    assert solution.flow.exit_mach < 1.0

    assert solution.flow.exit_pressure == pytest.approx(
        solution.back_pressure,
        rel=1.0e-10,
    )


def test_unchoked_mass_flow_is_below_critical_mass_flow():
    unchoked = solve_fanno_flow_for_back_pressure(
        p0=P0,
        T0=T0,
        back_pressure=8.0e5,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    critical = solve_critical_fanno_flow(
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert (
        unchoked.flow.inlet.mass_flow
        < critical.inlet.mass_flow
    )


def test_lower_back_pressure_does_not_change_choked_mass_flow():
    first = solve_fanno_flow_for_back_pressure(
        p0=P0,
        T0=T0,
        back_pressure=5.0e5,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    second = solve_fanno_flow_for_back_pressure(
        p0=P0,
        T0=T0,
        back_pressure=1.013e5,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert first.regime == "choked"
    assert second.regime == "choked"

    assert first.flow.inlet.mass_flow == pytest.approx(
        second.flow.inlet.mass_flow,
        rel=1.0e-12,
    )


@pytest.mark.parametrize(
    "back_pressure",
    [0.0, -1.0, P0, 1.1 * P0, math.inf, math.nan],
)
def test_rejects_invalid_back_pressure(back_pressure):
    with pytest.raises(ValueError):
        solve_fanno_flow_for_back_pressure(
            p0=P0,
            T0=T0,
            back_pressure=back_pressure,
            gamma=GAMMA,
            gas_constant=GAS_CONSTANT,
            diameter=DIAMETER,
            length=LENGTH,
            darcy_friction_factor=FRICTION_FACTOR,
        )