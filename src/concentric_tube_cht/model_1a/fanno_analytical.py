"""Analytical relations for constant-property Fanno flow."""

import math
from scipy.optimize import brentq
from dataclasses import dataclass

def fanno_length_parameter(M: float, gamma: float) -> float:
    """Return the dimensionless distance from Mach ``M`` to sonic conditions.

    The returned quantity is ``f_D * L_star / D``, where ``f_D`` is the
    Darcy friction factor.

    Parameters
    ----------
    M
        Local Mach number.
    gamma
        Specific-heat ratio.

    Returns
    -------
    float
        Dimensionless Fanno length parameter.

    Raises
    ------
    ValueError
        If the Mach number or specific-heat ratio is outside its physical
        domain.
    """
    if not math.isfinite(M) or M <= 0.0:
        raise ValueError("Mach number must be finite and greater than zero.")

    if not math.isfinite(gamma) or gamma <= 1.0:
        raise ValueError(
            "Specific-heat ratio must be finite and greater than one."
        )

    M_squared = M**2

    algebraic_term = (
        (1.0 - M_squared)
        / (gamma * M_squared)
    )

    logarithmic_argument = (
        (gamma + 1.0) * M_squared
        / (2.0 + (gamma - 1.0) * M_squared)
    )

    logarithmic_term = (
        (gamma + 1.0)
        / (2.0 * gamma)
        * math.log(logarithmic_argument)
    )

    return algebraic_term + logarithmic_term

def fanno_temperature_ratio(M: float, gamma: float) -> float:
    """Return the Fanno static-temperature ratio T / T_star."""
    if not math.isfinite(M) or M <= 0.0:
        raise ValueError("Mach number must be finite and greater than zero.")

    if not math.isfinite(gamma) or gamma <= 1.0:
        raise ValueError(
            "Specific-heat ratio must be finite and greater than one."
        )

    return (gamma + 1.0) / (
        2.0 + (gamma - 1.0) * M**2
    )


def fanno_pressure_ratio(M: float, gamma: float) -> float:
    """Return the Fanno static-pressure ratio p / p_star."""
    temperature_ratio = fanno_temperature_ratio(M, gamma)

    return math.sqrt(temperature_ratio) / M

def fanno_density_ratio(M: float, gamma: float) -> float:
    """Return the Fanno density ratio rho / rho_star."""
    temperature_ratio = fanno_temperature_ratio(M, gamma)

    return 1.0 / (M * math.sqrt(temperature_ratio))


def fanno_velocity_ratio(M: float, gamma: float) -> float:
    """Return the Fanno velocity ratio u / u_star."""
    temperature_ratio = fanno_temperature_ratio(M, gamma)

    return M * math.sqrt(temperature_ratio)


def fanno_stagnation_pressure_ratio(
    M: float,
    gamma: float,
) -> float:
    """Return the Fanno stagnation-pressure ratio p0 / p0_star."""
    temperature_ratio = fanno_temperature_ratio(M, gamma)

    exponent = (
        (gamma + 1.0)
        / (2.0 * (gamma - 1.0))
    )

    return (
        (1.0 / M)
        * (1.0 / temperature_ratio) ** exponent
    )

class SonicLimitError(RuntimeError):
    """Raised when a trial inlet state reaches Mach 1 before the exit."""


def solve_subsonic_exit_mach(
    M_in: float,
    gamma: float,
    darcy_friction_factor: float,
    length: float,
    diameter: float,
) -> float:
    """Return the subsonic exit Mach number for a prescribed inlet Mach.

    The function never solves beyond Mach 1.
    """
    values = (
        M_in,
        gamma,
        darcy_friction_factor,
        length,
        diameter,
    )

    if not all(math.isfinite(value) for value in values):
        raise ValueError("All inputs must be finite.")

    if not 0.0 < M_in < 1.0:
        raise ValueError("Inlet Mach number must satisfy 0 < M_in < 1.")

    if gamma <= 1.0:
        raise ValueError("Specific-heat ratio must be greater than one.")

    if darcy_friction_factor < 0.0 or length < 0.0:
        raise ValueError("Friction factor and length cannot be negative.")

    if diameter <= 0.0:
        raise ValueError("Diameter must be positive.")

    friction_length = (
        darcy_friction_factor * length / diameter
    )

    inlet_fanno_parameter = fanno_length_parameter(
        M=M_in,
        gamma=gamma,
    )

    remaining_fanno_parameter = (
        inlet_fanno_parameter - friction_length
    )

    tolerance = 1.0e-12

    if remaining_fanno_parameter < -tolerance:
        sonic_location = (
            diameter
            * inlet_fanno_parameter
            / darcy_friction_factor
        )

        raise SonicLimitError(
            "The trial inlet Mach number reaches M = 1 "
            f"at x = {sonic_location:.8g} m, before the tube exit."
        )

    if abs(remaining_fanno_parameter) <= tolerance:
        return 1.0

    if friction_length == 0.0:
        return M_in

    def residual(M_out: float) -> float:
        return (
            fanno_length_parameter(M_out, gamma)
            - remaining_fanno_parameter
        )

    return brentq(
        residual,
        M_in,
        1.0 - 1.0e-10,
        xtol=1.0e-12,
        rtol=1.0e-12,
    )

@dataclass(frozen=True)
class FannoStationRatios:
    """Property ratios between downstream and upstream Fanno states."""

    temperature: float
    pressure: float
    density: float
    velocity: float
    stagnation_pressure: float


def subsonic_fanno_station_ratios(
    M_upstream: float,
    M_downstream: float,
    gamma: float,
) -> FannoStationRatios:
    """Return downstream-to-upstream property ratios.

    Both states must lie on the same subsonic Fanno branch, with friction
    accelerating the flow toward Mach 1.
    """
    if not 0.0 < M_upstream <= M_downstream <= 1.0:
        raise ValueError(
            "Mach numbers must satisfy "
            "0 < M_upstream <= M_downstream <= 1."
        )

    temperature_ratio = (
        fanno_temperature_ratio(M_downstream, gamma)
        / fanno_temperature_ratio(M_upstream, gamma)
    )

    pressure_ratio = (
        fanno_pressure_ratio(M_downstream, gamma)
        / fanno_pressure_ratio(M_upstream, gamma)
    )

    density_ratio = (
        fanno_density_ratio(M_downstream, gamma)
        / fanno_density_ratio(M_upstream, gamma)
    )

    velocity_ratio = (
        fanno_velocity_ratio(M_downstream, gamma)
        / fanno_velocity_ratio(M_upstream, gamma)
    )

    stagnation_pressure_ratio = (
        fanno_stagnation_pressure_ratio(M_downstream, gamma)
        / fanno_stagnation_pressure_ratio(M_upstream, gamma)
    )

    return FannoStationRatios(
        temperature=temperature_ratio,
        pressure=pressure_ratio,
        density=density_ratio,
        velocity=velocity_ratio,
        stagnation_pressure=stagnation_pressure_ratio,
    )

def solve_critical_inlet_mach(
    gamma: float,
    darcy_friction_factor: float,
    length: float,
    diameter: float,
) -> float:
    """Return the inlet Mach number that produces M = 1 at the exit."""
    values = (
        gamma,
        darcy_friction_factor,
        length,
        diameter,
    )

    if not all(math.isfinite(value) for value in values):
        raise ValueError("All inputs must be finite.")

    if gamma <= 1.0:
        raise ValueError("Specific-heat ratio must be greater than one.")

    if darcy_friction_factor < 0.0 or length < 0.0:
        raise ValueError("Friction factor and length cannot be negative.")

    if diameter <= 0.0:
        raise ValueError("Diameter must be positive.")

    friction_length = (
        darcy_friction_factor * length / diameter
    )

    if friction_length == 0.0:
        return 1.0

    def residual(M_in: float) -> float:
        return (
            fanno_length_parameter(M_in, gamma)
            - friction_length
        )

    return brentq(
        residual,
        1.0e-8,
        1.0 - 1.0e-10,
        xtol=1.0e-12,
        rtol=1.0e-12,
    )