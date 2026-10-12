import math
import pytest

from concentric_tube_cht.model_1a.fanno_analytical import (
    FannoStationRatios,
    SonicLimitError,
    fanno_density_ratio,
    fanno_length_parameter,
    fanno_pressure_ratio,
    fanno_stagnation_pressure_ratio,
    fanno_temperature_ratio,
    fanno_velocity_ratio,
    solve_subsonic_exit_mach,
    subsonic_fanno_station_ratios,
)

def test_sonic_state_has_zero_remaining_length():
    result = fanno_length_parameter(M=1.0, gamma=1.4)

    assert result == pytest.approx(0.0, abs=1.0e-14)


@pytest.mark.parametrize(
    ("M", "expected"),
    [
        (0.5, 1.0690603127182559),
        (2.0, 0.3049965025814798),
    ],
)
def test_matches_reference_values(M, expected):
    result = fanno_length_parameter(M=M, gamma=1.4)

    assert result == pytest.approx(expected, rel=1.0e-12)


@pytest.mark.parametrize("M", [0.5, 0.8, 1.2, 2.0])
def test_length_parameter_is_positive_away_from_sonic_state(M):
    result = fanno_length_parameter(M=M, gamma=1.4)

    assert result > 0.0


def test_converges_to_zero_from_both_sides_of_sonic_state():
    subsonic_result = fanno_length_parameter(M=0.999, gamma=1.4)
    supersonic_result = fanno_length_parameter(M=1.001, gamma=1.4)

    assert subsonic_result > 0.0
    assert supersonic_result > 0.0
    assert subsonic_result < 1.0e-5
    assert supersonic_result < 1.0e-5


@pytest.mark.parametrize("M", [0.0, -0.5, math.inf, math.nan])
def test_rejects_invalid_mach_numbers(M):
    with pytest.raises(ValueError):
        fanno_length_parameter(M=M, gamma=1.4)


@pytest.mark.parametrize("gamma", [1.0, 0.9, math.inf, math.nan])
def test_rejects_invalid_specific_heat_ratios(gamma):
    with pytest.raises(ValueError):
        fanno_length_parameter(M=0.5, gamma=gamma)


@pytest.mark.parametrize(
    ("M", "expected_temperature_ratio", "expected_pressure_ratio"),
    [
        (0.5, 1.1428571428571428, 2.138089935299395),
        (1.0, 1.0, 1.0),
        (2.0, 0.6666666666666666, 0.408248290463863),
    ],
)
def test_fanno_temperature_and_pressure_ratios(
    M,
    expected_temperature_ratio,
    expected_pressure_ratio,
):
    temperature_ratio = fanno_temperature_ratio(M, gamma=1.4)
    pressure_ratio = fanno_pressure_ratio(M, gamma=1.4)

    assert temperature_ratio == pytest.approx(
        expected_temperature_ratio,
        rel=1.0e-12,
    )

    assert pressure_ratio == pytest.approx(
        expected_pressure_ratio,
        rel=1.0e-12,
    )

@pytest.mark.parametrize(
    (
        "M",
        "expected_density_ratio",
        "expected_velocity_ratio",
        "expected_stagnation_pressure_ratio",
    ),
    [
        (
            0.5,
            1.8708286933869707,
            0.5345224838248488,
            1.33984375,
        ),
        (
            1.0,
            1.0,
            1.0,
            1.0,
        ),
        (
            2.0,
            0.6123724356957945,
            1.632993161855452,
            1.6875,
        ),
    ],
)
def test_remaining_fanno_property_ratios(
    M,
    expected_density_ratio,
    expected_velocity_ratio,
    expected_stagnation_pressure_ratio,
):
    density_ratio = fanno_density_ratio(M, gamma=1.4)
    velocity_ratio = fanno_velocity_ratio(M, gamma=1.4)
    stagnation_pressure_ratio = (
        fanno_stagnation_pressure_ratio(M, gamma=1.4)
    )

    assert density_ratio == pytest.approx(
        expected_density_ratio,
        rel=1.0e-12,
    )

    assert velocity_ratio == pytest.approx(
        expected_velocity_ratio,
        rel=1.0e-12,
    )

    assert stagnation_pressure_ratio == pytest.approx(
        expected_stagnation_pressure_ratio,
        rel=1.0e-12,
    )


@pytest.mark.parametrize("M", [0.5, 0.8, 1.0, 1.2, 2.0])
def test_density_velocity_ratios_satisfy_continuity(M):
    density_ratio = fanno_density_ratio(M, gamma=1.4)
    velocity_ratio = fanno_velocity_ratio(M, gamma=1.4)

    assert density_ratio * velocity_ratio == pytest.approx(
        1.0,
        rel=1.0e-12,
    )

def test_subsonic_exit_mach_satisfies_fanno_length_equation():
    M_in = 0.5
    gamma = 1.4
    friction_factor = 0.02
    length = 0.020
    diameter = 0.010

    M_out = solve_subsonic_exit_mach(
        M_in=M_in,
        gamma=gamma,
        darcy_friction_factor=friction_factor,
        length=length,
        diameter=diameter,
    )

    residual = (
        fanno_length_parameter(M_in, gamma)
        - fanno_length_parameter(M_out, gamma)
        - friction_factor * length / diameter
    )

    assert M_out == pytest.approx(
        0.5049924572060341,
        rel=1.0e-10,
    )
    assert residual == pytest.approx(0.0, abs=1.0e-12)
    assert M_in < M_out < 1.0


def test_zero_friction_preserves_mach_number():
    M_out = solve_subsonic_exit_mach(
        M_in=0.5,
        gamma=1.4,
        darcy_friction_factor=0.0,
        length=0.020,
        diameter=0.010,
    )

    assert M_out == pytest.approx(0.5)


def test_returns_sonic_mach_when_sonic_state_is_at_exit():
    M_in = 0.5
    gamma = 1.4
    friction_factor = 0.02
    diameter = 0.010

    critical_length = (
        diameter
        * fanno_length_parameter(M_in, gamma)
        / friction_factor
    )

    M_out = solve_subsonic_exit_mach(
        M_in=M_in,
        gamma=gamma,
        darcy_friction_factor=friction_factor,
        length=critical_length,
        diameter=diameter,
    )

    assert M_out == pytest.approx(1.0)


def test_rejects_trial_state_that_chokes_before_exit():
    M_in = 0.5
    gamma = 1.4
    friction_factor = 0.02
    diameter = 0.010

    critical_length = (
        diameter
        * fanno_length_parameter(M_in, gamma)
        / friction_factor
    )

    with pytest.raises(SonicLimitError):
        solve_subsonic_exit_mach(
            M_in=M_in,
            gamma=gamma,
            darcy_friction_factor=friction_factor,
            length=1.01 * critical_length,
            diameter=diameter,
        )

def test_subsonic_station_ratios_match_reference_values():
    ratios = subsonic_fanno_station_ratios(
        M_upstream=0.5,
        M_downstream=0.5049924572060341,
        gamma=1.4,
    )

    assert ratios.temperature == pytest.approx(
        0.9990452206966373,
        rel=1.0e-12,
    )
    assert ratios.pressure == pytest.approx(
        0.9896410155049582,
        rel=1.0e-12,
    )
    assert ratios.density == pytest.approx(
        0.9905868072867396,
        rel=1.0e-12,
    )
    assert ratios.velocity == pytest.approx(
        1.009502642922374,
        rel=1.0e-12,
    )
    assert ratios.stagnation_pressure == pytest.approx(
        0.9929552431298154,
        rel=1.0e-12,
    )


def test_station_ratios_satisfy_mass_conservation():
    ratios = subsonic_fanno_station_ratios(
        M_upstream=0.5,
        M_downstream=0.7,
        gamma=1.4,
    )

    assert ratios.density * ratios.velocity == pytest.approx(
        1.0,
        rel=1.0e-12,
    )


def test_station_ratios_satisfy_perfect_gas_relation():
    ratios = subsonic_fanno_station_ratios(
        M_upstream=0.5,
        M_downstream=0.7,
        gamma=1.4,
    )

    assert ratios.pressure == pytest.approx(
        ratios.density * ratios.temperature,
        rel=1.0e-12,
    )


def test_identical_mach_numbers_produce_unit_ratios():
    ratios = subsonic_fanno_station_ratios(
        M_upstream=0.5,
        M_downstream=0.5,
        gamma=1.4,
    )

    assert ratios.temperature == pytest.approx(1.0)
    assert ratios.pressure == pytest.approx(1.0)
    assert ratios.density == pytest.approx(1.0)
    assert ratios.velocity == pytest.approx(1.0)
    assert ratios.stagnation_pressure == pytest.approx(1.0)