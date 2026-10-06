# Problem Definition

## Objective

This project develops a hierarchy of analytical and numerical models for conjugate heat transfer between a hot gas and a counter-flowing liquid coolant in concentric copper tubes.

The initial model is intended as a transparent engineering screening calculation. It is not a validated design prediction.

## Geometry

The system consists of two concentric copper tubes with a common heated length of 20 mm.

| Parameter | Value |
|---|---:|
| Heated length | 20 mm |
| Inner-tube inner diameter | 10 mm |
| Inner-tube outer diameter | 12 mm |
| Outer-tube inner diameter | 14 mm |
| Outer-tube outer diameter | 16 mm |

The hot gas flows through the inner tube. Liquid water flows through the annular passage between the two tubes in the opposite direction.

## Operating Conditions

### Hot-gas stream

| Parameter | Value |
|---|---:|
| Composition | 50 mol% H₂O, 50 mol% CO₂ |
| Inlet temperature | 1000 K |
| Inlet pressure | 10 bar |
| Downstream pressure | Atmospheric |
| Flow direction | Inner-tube stream |

Gas-mixture properties are evaluated using Cantera with the GRI-Mech 3.0 mechanism.

### Coolant stream

| Parameter | Value |
|---|---:|
| Fluid | Liquid water |
| Inlet temperature | 300 K |
| Inlet pressure | 15 bar |
| Mass-flow rate | 0.5 kg/s |
| Flow direction | Annular counter-flow |

Water properties are evaluated using CoolProp.

## Model Hierarchy

### Model 0 — Constant-property screening model

Model 0 is a steady, lumped thermal-resistance model using constant thermophysical properties. It estimates:

- Flow areas and hydraulic diameter
- Coolant pressure drop
- Ideal gas choking and mass-flow screening
- Gas-side and coolant-side heat-transfer coefficients
- Radial thermal resistances
- Heat-transfer rate and wall heat fluxes
- Interface temperatures
- Approximate outlet temperatures
- Coolant saturation margin

The ideal choked gas mass flow is used provisionally in the gas-side heat-transfer calculation.

### Model 1 — Coupled one-dimensional model

The planned next stage will introduce a coupled axial solution with variable properties and compressible gas-flow effects. Model 1 is outside the current repository milestone and has not yet been implemented.

## Scope Statement

Model 0 provides a reproducible baseline for checking orders of magnitude, identifying dominant thermal resistances, and guiding subsequent model development. Its results should not be interpreted as final hardware performance or validated operating limits.