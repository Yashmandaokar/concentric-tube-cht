"""Reservoir-to-tube inlet relations for a constant perfect gas."""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class InletState:
    """Thermodynamic and flow state at the tube inlet."""

    mach: float
    pressure: float
    temperature: float
    density: float
    speed_of_sound: float
    velocity: float
    mass_flow: float


def reservoir_to_inlet_state(
    M_in: float,
    p0: float,
    T0: float,
    gamma: float,
    gas_constant: float,
    area: float,
) -> InletState:
    """Calculate the tube-inlet state from reservoir stagnation conditions."""
    values = (M_in, p0, T0, gamma, gas_constant, area)

    if not all(math.isfinite(value) for value in values):
        raise ValueError("All inputs must be finite.")

    if not 0.0 < M_in < 1.0:
        raise ValueError("The trial inlet Mach number must satisfy 0 < M < 1.")

    if p0 <= 0.0 or T0 <= 0.0:
        raise ValueError("Reservoir pressure and temperature must be positive.")

    if gamma <= 1.0:
        raise ValueError("Specific-heat ratio must be greater than one.")

    if gas_constant <= 0.0 or area <= 0.0:
        raise ValueError("Gas constant and flow area must be positive.")

    stagnation_factor = (
        1.0
        + 0.5 * (gamma - 1.0) * M_in**2
    )

    temperature = T0 / stagnation_factor

    pressure = p0 / (
        stagnation_factor ** (gamma / (gamma - 1.0))
    )

    density = pressure / (gas_constant * temperature)

    speed_of_sound = math.sqrt(
        gamma * gas_constant * temperature
    )

    velocity = M_in * speed_of_sound
    mass_flow = density * velocity * area

    return InletState(
        mach=M_in,
        pressure=pressure,
        temperature=temperature,
        density=density,
        speed_of_sound=speed_of_sound,
        velocity=velocity,
        mass_flow=mass_flow,
    )