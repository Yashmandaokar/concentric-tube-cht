"""Mach-based analytical reference solver for Model 1A."""

import math
from dataclasses import dataclass
from typing import Literal
from scipy.optimize import brentq

from .fanno_analytical import (
    solve_critical_inlet_mach,
    solve_subsonic_exit_mach,
    subsonic_fanno_station_ratios,
)
from .reservoir_relations import (
    InletState,
    reservoir_to_inlet_state,
)


@dataclass(frozen=True)
class AnalyticalFannoSolution:
    """Analytical Fanno solution for a prescribed inlet Mach number."""

    inlet: InletState
    exit_mach: float
    exit_pressure: float
    exit_temperature: float
    exit_density: float
    exit_speed_of_sound: float
    exit_velocity: float
    exit_stagnation_pressure: float
    choked_at_exit: bool

@dataclass(frozen=True)
class BackPressureSolution:
    """Physical Fanno solution selected using the imposed back pressure."""

    regime: Literal["unchoked", "choked"]
    flow: AnalyticalFannoSolution
    back_pressure: float
    critical_back_pressure: float

def solve_trial_fanno_flow(
    M_in: float,
    p0: float,
    T0: float,
    gamma: float,
    gas_constant: float,
    diameter: float,
    length: float,
    darcy_friction_factor: float,
) -> AnalyticalFannoSolution:
    """Solve the tube flow for a prescribed subsonic inlet Mach number."""
    if diameter <= 0.0:
        raise ValueError("Diameter must be positive.")

    area = math.pi * diameter**2 / 4.0

    inlet = reservoir_to_inlet_state(
        M_in=M_in,
        p0=p0,
        T0=T0,
        gamma=gamma,
        gas_constant=gas_constant,
        area=area,
    )

    exit_mach = solve_subsonic_exit_mach(
        M_in=M_in,
        gamma=gamma,
        darcy_friction_factor=darcy_friction_factor,
        length=length,
        diameter=diameter,
    )

    ratios = subsonic_fanno_station_ratios(
        M_upstream=M_in,
        M_downstream=exit_mach,
        gamma=gamma,
    )

    exit_temperature = (
        inlet.temperature * ratios.temperature
    )

    exit_pressure = (
        inlet.pressure * ratios.pressure
    )

    exit_density = (
        inlet.density * ratios.density
    )

    exit_velocity = (
        inlet.velocity * ratios.velocity
    )

    exit_speed_of_sound = math.sqrt(
        gamma * gas_constant * exit_temperature
    )

    exit_stagnation_pressure = (
        p0 * ratios.stagnation_pressure
    )

    choked_at_exit = math.isclose(
        exit_mach,
        1.0,
        rel_tol=0.0,
        abs_tol=1.0e-10,
    )

    return AnalyticalFannoSolution(
        inlet=inlet,
        exit_mach=exit_mach,
        exit_pressure=exit_pressure,
        exit_temperature=exit_temperature,
        exit_density=exit_density,
        exit_speed_of_sound=exit_speed_of_sound,
        exit_velocity=exit_velocity,
        exit_stagnation_pressure=exit_stagnation_pressure,
        choked_at_exit=choked_at_exit,
    )

def solve_critical_fanno_flow(
    p0: float,
    T0: float,
    gamma: float,
    gas_constant: float,
    diameter: float,
    length: float,
    darcy_friction_factor: float,
) -> AnalyticalFannoSolution:
    """Return the critical solution with Mach 1 at the tube exit."""
    critical_inlet_mach = solve_critical_inlet_mach(
        gamma=gamma,
        darcy_friction_factor=darcy_friction_factor,
        length=length,
        diameter=diameter,
    )

    return solve_trial_fanno_flow(
        M_in=critical_inlet_mach,
        p0=p0,
        T0=T0,
        gamma=gamma,
        gas_constant=gas_constant,
        diameter=diameter,
        length=length,
        darcy_friction_factor=darcy_friction_factor,
    )

def solve_fanno_flow_for_back_pressure(
    p0: float,
    T0: float,
    back_pressure: float,
    gamma: float,
    gas_constant: float,
    diameter: float,
    length: float,
    darcy_friction_factor: float,
) -> BackPressureSolution:
    """Select and solve the choked or unchoked Fanno-flow branch."""
    if not math.isfinite(back_pressure) or back_pressure <= 0.0:
        raise ValueError("Back pressure must be finite and positive.")

    if back_pressure >= p0:
        raise ValueError(
            "Back pressure must be below reservoir pressure "
            "for forward flow."
        )

    critical_flow = solve_critical_fanno_flow(
        p0=p0,
        T0=T0,
        gamma=gamma,
        gas_constant=gas_constant,
        diameter=diameter,
        length=length,
        darcy_friction_factor=darcy_friction_factor,
    )

    critical_back_pressure = critical_flow.exit_pressure

    if back_pressure <= critical_back_pressure:
        return BackPressureSolution(
            regime="choked",
            flow=critical_flow,
            back_pressure=back_pressure,
            critical_back_pressure=critical_back_pressure,
        )

    critical_inlet_mach = critical_flow.inlet.mach

    def exit_pressure_residual(M_in: float) -> float:
        trial_flow = solve_trial_fanno_flow(
            M_in=M_in,
            p0=p0,
            T0=T0,
            gamma=gamma,
            gas_constant=gas_constant,
            diameter=diameter,
            length=length,
            darcy_friction_factor=darcy_friction_factor,
        )

        return trial_flow.exit_pressure - back_pressure

    unchoked_inlet_mach = brentq(
        exit_pressure_residual,
        1.0e-8,
        critical_inlet_mach,
        xtol=1.0e-12,
        rtol=1.0e-12,
    )

    unchoked_flow = solve_trial_fanno_flow(
        M_in=unchoked_inlet_mach,
        p0=p0,
        T0=T0,
        gamma=gamma,
        gas_constant=gas_constant,
        diameter=diameter,
        length=length,
        darcy_friction_factor=darcy_friction_factor,
    )

    return BackPressureSolution(
        regime="unchoked",
        flow=unchoked_flow,
        back_pressure=back_pressure,
        critical_back_pressure=critical_back_pressure,
    )