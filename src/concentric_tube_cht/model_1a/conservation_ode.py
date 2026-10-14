"""Conservation-equation ODE formulation for Model 1A."""

import math
from dataclasses import dataclass
from scipy.integrate import solve_ivp
import numpy as np


class SonicPointError(RuntimeError):
    """Raised when the ODE evaluation reaches the sonic singularity."""


@dataclass(frozen=True)
class FannoODEParameters:
    """Constant parameters required by the p-T ODE system."""

    mass_flow: float
    area: float
    gamma: float
    gas_constant: float
    diameter: float
    darcy_friction_factor: float

@dataclass(frozen=True)
class FannoODEResult:
    """Axial solution produced by the p-T ODE solver."""

    x: np.ndarray
    pressure: np.ndarray
    temperature: np.ndarray
    density: np.ndarray
    velocity: np.ndarray
    mach: np.ndarray
    terminated_by_sonic_event: bool

def fanno_pt_rhs(
    x: float,
    state: np.ndarray,
    parameters: FannoODEParameters,
) -> np.ndarray:
    """Return dp/dx and dT/dx for adiabatic constant-area flow."""
    pressure, temperature = state

    if not math.isfinite(pressure) or pressure <= 0.0:
        raise ValueError("Pressure must be finite and positive.")

    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("Temperature must be finite and positive.")

    density = (
        pressure
        / (parameters.gas_constant * temperature)
    )

    velocity = (
        parameters.mass_flow
        / (density * parameters.area)
    )

    speed_of_sound = math.sqrt(
        parameters.gamma
        * parameters.gas_constant
        * temperature
    )

    mach = velocity / speed_of_sound

    if mach >= 1.0 - 1.0e-10:
        raise SonicPointError(
            "The p-T ODE system cannot be evaluated at M = 1."
        )

    cp = (
        parameters.gamma
        * parameters.gas_constant
        / (parameters.gamma - 1.0)
    )

    coefficient_matrix = np.array(
        [
            [
                -(velocity**2) / pressure,
                cp + velocity**2 / temperature,
            ],
            [
                1.0 - density * velocity**2 / pressure,
                density * velocity**2 / temperature,
            ],
        ],
        dtype=float,
    )

    source_vector = np.array(
        [
            0.0,
            (
                -parameters.darcy_friction_factor
                * density
                * velocity**2
                / (2.0 * parameters.diameter)
            ),
        ],
        dtype=float,
    )

    dp_dx, dT_dx = np.linalg.solve(
        coefficient_matrix,
        source_vector,
    )

    return np.array([dp_dx, dT_dx])

def integrate_fanno_pt(
    inlet_pressure: float,
    inlet_temperature: float,
    length: float,
    parameters: FannoODEParameters,
    relative_tolerance: float = 1.0e-10,
    sonic_mach_limit: float = 0.999,
) -> FannoODEResult:
    """Integrate the p-T conservation equations through the tube."""
    if length <= 0.0 or not math.isfinite(length):
        raise ValueError("Tube length must be finite and positive.")

    if not 0.0 < sonic_mach_limit < 1.0:
        raise ValueError(
            "Sonic Mach limit must satisfy 0 < limit < 1."
        )

    def rhs(x: float, state: np.ndarray) -> np.ndarray:
        return fanno_pt_rhs(
            x=x,
            state=state,
            parameters=parameters,
        )

    def sonic_event(
        x: float,
        state: np.ndarray,
    ) -> float:
        pressure, temperature = state

        density = (
            pressure
            / (parameters.gas_constant * temperature)
        )

        velocity = (
            parameters.mass_flow
            / (density * parameters.area)
        )

        speed_of_sound = math.sqrt(
            parameters.gamma
            * parameters.gas_constant
            * temperature
        )

        mach = velocity / speed_of_sound

        return sonic_mach_limit - mach

    sonic_event.terminal = True
    sonic_event.direction = -1.0

    integration = solve_ivp(
        fun=rhs,
        t_span=(0.0, length),
        y0=np.array([
            inlet_pressure,
            inlet_temperature,
        ]),
        method="DOP853",
        rtol=relative_tolerance,
        atol=np.array([1.0e-5, 1.0e-9]),
        max_step=length / 100.0,
        events=sonic_event,
    )

    if not integration.success:
        raise RuntimeError(
            f"ODE integration failed: {integration.message}"
        )

    x = integration.t
    pressure = integration.y[0]
    temperature = integration.y[1]

    density = (
        pressure
        / (parameters.gas_constant * temperature)
    )

    velocity = (
        parameters.mass_flow
        / (density * parameters.area)
    )

    speed_of_sound = np.sqrt(
        parameters.gamma
        * parameters.gas_constant
        * temperature
    )

    mach = velocity / speed_of_sound

    terminated_by_sonic_event = (
        len(integration.t_events[0]) > 0
    )

    return FannoODEResult(
        x=x,
        pressure=pressure,
        temperature=temperature,
        density=density,
        velocity=velocity,
        mach=mach,
        terminated_by_sonic_event=terminated_by_sonic_event,
    )