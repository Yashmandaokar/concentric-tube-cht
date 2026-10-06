# Concentric-Tube Conjugate Heat-Transfer Model

A reproducible engineering portfolio project investigating conjugate heat transfer between a hot H₂O–CO₂ gas mixture and counter-flowing liquid water in concentric copper tubes.

The project follows a staged model hierarchy, beginning with a transparent screening calculation and progressing toward a coupled one-dimensional compressible-flow model.

> **Current status:** Model 0 is complete and frozen. It is a screening-level baseline, not a validated design prediction.

## Engineering Problem

The system consists of two concentric copper tubes with a heated length of 20 mm.

| Parameter | Value |
|---|---:|
| Inner-tube ID / OD | 10 mm / 12 mm |
| Outer-tube ID / OD | 14 mm / 16 mm |
| Hot gas | 50 mol% H₂O, 50 mol% CO₂ |
| Gas inlet condition | 1000 K, 10 bar |
| Gas downstream pressure | Atmospheric |
| Coolant | Liquid water |
| Coolant inlet condition | 300 K, 15 bar |
| Coolant mass-flow rate | 0.5 kg/s |
| Flow arrangement | Counter-flow |

Gas properties are evaluated with Cantera using GRI-Mech 3.0. Water properties are evaluated with CoolProp.

## Model Hierarchy

### Model 0 — Constant-property baseline

The completed Model 0 notebook includes:

- Flow geometry and characteristic dimensions
- Gas and coolant thermophysical properties
- Coolant Reynolds number, friction factor, and pressure drop
- Ideal gas choking screening
- Provisional ideal choked gas mass flow
- Gas-side and coolant-side Gnielinski correlations
- Cylindrical copper conduction
- Radial thermal-resistance network
- Heat rate, wall heat fluxes, and interface temperatures
- Approximate fluid outlet temperatures
- Energy-conservation assertions
- Coolant saturation-temperature screening

### Model 1 — Planned

The next model level will introduce an axially coupled one-dimensional solution with variable properties and compressible gas-flow effects. It has not yet been implemented.

## Principal Model 0 Results

| Quantity | Result |
|---|---:|
| Ideal gas mass-flow estimate | 0.09868 kg/s |
| Reference gas Mach number | 0.591 |
| Gas-side HTC | 3924 W/(m²·K) |
| Coolant-side HTC | 57,647 W/(m²·K) |
| Total thermal resistance | 0.43241 K/W |
| Gas-side resistance contribution | 93.80% |
| Heat-transfer rate | 1.619 kW |
| Gas-side wall temperature | 343.43 K |
| Coolant-side wall temperature | 337.24 K |
| Approximate gas outlet temperature | 989.36 K |
| Approximate coolant outlet temperature | 300.78 K |

At the assumed coolant pressure of 15 bar, the calculated coolant-side wall temperature remains below the water saturation temperature of approximately 471 K.

## Engineering Interpretation

Model 0 identifies gas-side convection as the dominant thermal resistance. However, the gas passage has \(L/D=2\) and a reference Mach number of approximately 0.59. The assumptions behind a fully developed incompressible heat-transfer correlation are therefore weak on the gas side.

The ideal choking calculation is used only to obtain a provisional mass-flow estimate. It does not predict flow through the cooled, frictional, constant-area tube.

For these reasons, the results establish orders of magnitude and guide further model development; they should not be interpreted as final hardware performance.

## Installation

Python 3.11 or a compatible Python version is recommended.

Clone the repository and enter its directory:

```bash
git clone https://github.com/<your-username>/concentric-tube-cht.git
cd concentric-tube-cht
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Windows:

```powershell
.venv\Scripts\Activate.ps1
```

On macOS or Linux:

```bash
source .venv/bin/activate
```

Install the dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Launch Jupyter:

```bash
jupyter notebook
```

Then open:

```text
notebooks/model_0/model_0_constant_property_baseline.ipynb
```

Run the notebook from top to bottom to reproduce the Model 0 calculations.

## Documentation

- [Problem definition](docs/problem_definition.md)
- [Assumptions and limitations](docs/assumptions_and_limitations.md)
- [Model 0 results summary](results/model_0/summary.md)

## Disclaimer

This repository is an independent engineering portfolio project. It is intended for educational and model-development purposes and does not represent certified hardware analysis or validated design data.