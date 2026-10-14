import math

import numpy as np
import pytest

from concentric_tube_cht.model_1a.conservation_ode import (
    FannoODEParameters,
    SonicPointError,
    fanno_pt_rhs,
    integrate_fanno_pt,
)
from concentric_tube_cht.model_1a.reservoir_relations import (
    reservoir_to_inlet_state,
)

from concentric_tube_cht.model_1a.analytical_solver import (
    solve_critical_fanno_flow,
    solve_fanno_flow_for_back_pressure,
    solve_trial_fanno_flow,
)

from concentric_tube_cht.model_1a.fanno_analytical import (
    fanno_length_parameter,
    subsonic_fanno_station_ratios,
)

P0 = 10.0e5
T0 = 1000.0
GAMMA = 1.2105
GAS_CONSTANT = 268.10
DIAMETER = 0.010
FRICTION_FACTOR = 0.01430
LENGTH = 0.020
AREA = math.pi * DIAMETER**2 / 4.0


def make_inlet_state():
    return reservoir_to_inlet_state(
        M_in=0.5,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        area=AREA,
    )


def make_parameters(mass_flow):
    return FannoODEParameters(
        mass_flow=mass_flow,
        area=AREA,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        darcy_friction_factor=FRICTION_FACTOR,
    )


def test_rhs_satisfies_energy_and_momentum_equations():
    inlet = make_inlet_state()
    parameters = make_parameters(inlet.mass_flow)

    dp_dx, dT_dx = fanno_pt_rhs(
        x=0.0,
        state=np.array([
            inlet.pressure,
            inlet.temperature,
        ]),
        parameters=parameters,
    )

    density = inlet.density
    velocity = inlet.velocity

    du_dx = velocity * (
        dT_dx / inlet.temperature
        - dp_dx / inlet.pressure
    )

    cp = (
        GAMMA
        * GAS_CONSTANT
        / (GAMMA - 1.0)
    )

    energy_residual = cp * dT_dx + velocity * du_dx

    momentum_residual = (
        dp_dx
        + density * velocity * du_dx
        + FRICTION_FACTOR
        * density
        * velocity**2
        / (2.0 * DIAMETER)
    )

    energy_scale = (
        abs(cp * dT_dx)
        + abs(velocity * du_dx)
    )

    momentum_scale = (
        abs(dp_dx)
        + abs(density * velocity * du_dx)
        + abs(
            FRICTION_FACTOR
            * density
            * velocity**2
            / (2.0 * DIAMETER)
        )
    )

    assert abs(energy_residual) / energy_scale < 1.0e-12
    assert abs(momentum_residual) / momentum_scale < 1.0e-12


def test_subsonic_friction_has_expected_derivative_signs():
    inlet = make_inlet_state()
    parameters = make_parameters(inlet.mass_flow)

    dp_dx, dT_dx = fanno_pt_rhs(
        x=0.0,
        state=np.array([
            inlet.pressure,
            inlet.temperature,
        ]),
        parameters=parameters,
    )

    du_dx = inlet.velocity * (
        dT_dx / inlet.temperature
        - dp_dx / inlet.pressure
    )

    assert dp_dx < 0.0
    assert dT_dx < 0.0
    assert du_dx > 0.0


def test_rhs_rejects_sonic_state():
    temperature = 900.0
    pressure = 5.0e5

    density = pressure / (GAS_CONSTANT * temperature)
    speed_of_sound = math.sqrt(
        GAMMA * GAS_CONSTANT * temperature
    )

    sonic_mass_flow = density * speed_of_sound * AREA
    parameters = make_parameters(sonic_mass_flow)

    with pytest.raises(SonicPointError):
        fanno_pt_rhs(
            x=0.0,
            state=np.array([pressure, temperature]),
            parameters=parameters,
        )

def test_unchoked_integration_matches_analytical_exit_state():
    inlet = make_inlet_state()
    parameters = make_parameters(inlet.mass_flow)

    numerical = integrate_fanno_pt(
        inlet_pressure=inlet.pressure,
        inlet_temperature=inlet.temperature,
        length=LENGTH,
        parameters=parameters,
    )

    analytical = solve_trial_fanno_flow(
        M_in=0.5,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert not numerical.terminated_by_sonic_event

    assert numerical.x[-1] == pytest.approx(
        LENGTH,
        abs=1.0e-14,
    )

    assert numerical.mach[-1] == pytest.approx(
        analytical.exit_mach,
        rel=1.0e-9,
    )

    assert numerical.pressure[-1] == pytest.approx(
        analytical.exit_pressure,
        rel=1.0e-9,
    )

    assert numerical.temperature[-1] == pytest.approx(
        analytical.exit_temperature,
        rel=1.0e-9,
    )

def test_sonic_event_location_matches_analytical_fanno_result():
    M_in = 0.9
    sonic_mach_limit = 0.999

    inlet = reservoir_to_inlet_state(
        M_in=M_in,
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        area=AREA,
    )

    parameters = make_parameters(inlet.mass_flow)

    numerical = integrate_fanno_pt(
        inlet_pressure=inlet.pressure,
        inlet_temperature=inlet.temperature,
        length=LENGTH,
        parameters=parameters,
        sonic_mach_limit=sonic_mach_limit,
    )

    expected_event_location = (
        DIAMETER
        / FRICTION_FACTOR
        * (
            fanno_length_parameter(M_in, GAMMA)
            - fanno_length_parameter(
                sonic_mach_limit,
                GAMMA,
            )
        )
    )

    assert numerical.terminated_by_sonic_event

    assert numerical.x[-1] < LENGTH

    assert numerical.mach[-1] == pytest.approx(
        sonic_mach_limit,
        abs=1.0e-9,
    )

    assert numerical.x[-1] == pytest.approx(
        expected_event_location,
        rel=1.0e-7,
    )

def test_critical_mass_flow_approaches_sonic_state_at_exit():
    sonic_mach_limit = 0.999

    critical = solve_critical_fanno_flow(
        p0=P0,
        T0=T0,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    parameters = make_parameters(
        critical.inlet.mass_flow
    )

    numerical = integrate_fanno_pt(
        inlet_pressure=critical.inlet.pressure,
        inlet_temperature=critical.inlet.temperature,
        length=LENGTH,
        parameters=parameters,
        sonic_mach_limit=sonic_mach_limit,
    )

    expected_event_location = (
        DIAMETER
        / FRICTION_FACTOR
        * (
            fanno_length_parameter(
                critical.inlet.mach,
                GAMMA,
            )
            - fanno_length_parameter(
                sonic_mach_limit,
                GAMMA,
            )
        )
    )

    expected_ratios = subsonic_fanno_station_ratios(
        M_upstream=critical.inlet.mach,
        M_downstream=sonic_mach_limit,
        gamma=GAMMA,
    )

    expected_pressure = (
        critical.inlet.pressure
        * expected_ratios.pressure
    )

    expected_temperature = (
        critical.inlet.temperature
        * expected_ratios.temperature
    )

    assert numerical.terminated_by_sonic_event

    assert numerical.mach[-1] == pytest.approx(
        sonic_mach_limit,
        abs=1.0e-9,
    )

    assert numerical.x[-1] == pytest.approx(
        expected_event_location,
        rel=1.0e-7,
    )

    assert numerical.pressure[-1] == pytest.approx(
        expected_pressure,
        rel=1.0e-8,
    )

    assert numerical.temperature[-1] == pytest.approx(
        expected_temperature,
        rel=1.0e-8,
    )

    assert numerical.x[-1] < LENGTH

    assert LENGTH - numerical.x[-1] < 2.0e-6

def test_ode_matches_physical_unchoked_back_pressure_solution():
    back_pressure = 8.0e5

    reference = solve_fanno_flow_for_back_pressure(
        p0=P0,
        T0=T0,
        back_pressure=back_pressure,
        gamma=GAMMA,
        gas_constant=GAS_CONSTANT,
        diameter=DIAMETER,
        length=LENGTH,
        darcy_friction_factor=FRICTION_FACTOR,
    )

    assert reference.regime == "unchoked"

    parameters = make_parameters(
        reference.flow.inlet.mass_flow
    )

    numerical = integrate_fanno_pt(
        inlet_pressure=reference.flow.inlet.pressure,
        inlet_temperature=reference.flow.inlet.temperature,
        length=LENGTH,
        parameters=parameters,
    )

    assert not numerical.terminated_by_sonic_event

    assert numerical.x[-1] == pytest.approx(
        LENGTH,
        abs=1.0e-14,
    )

    assert numerical.pressure[-1] == pytest.approx(
        back_pressure,
        rel=1.0e-9,
    )

    assert numerical.pressure[-1] == pytest.approx(
        reference.flow.exit_pressure,
        rel=1.0e-9,
    )

    assert numerical.temperature[-1] == pytest.approx(
        reference.flow.exit_temperature,
        rel=1.0e-9,
    )

    assert numerical.mach[-1] == pytest.approx(
        reference.flow.exit_mach,
        rel=1.0e-9,
    )