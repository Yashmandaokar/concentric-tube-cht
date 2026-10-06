# Model 0 Results Summary

## Model Description

Model 0 is a steady, lumped, constant-property thermal-resistance model. It provides a screening-level baseline and should not be interpreted as a validated design prediction.

## Gas-Flow Screening

| Quantity | Result |
|---|---:|
| Ideal choked mass-flow rate | 0.09868 kg/s |
| Reference gas velocity | 336.9 m/s |
| Reference Mach number | 0.591 |
| Reference Reynolds number | 3.153 × 10⁵ |
| Gas-side heat-transfer coefficient | 3924 W/(m²·K) |

The gas mass flow is based on ideal isentropic choking. It is used provisionally for the Model 0 heat-transfer calculation and is not a coupled frictional-flow prediction.

## Coolant-Side Results

| Quantity | Result |
|---|---:|
| Coolant mass-flow rate | 0.5 kg/s |
| Coolant Reynolds number | 2.868 × 10⁴ |
| Coolant Nusselt number | 188.9 |
| Coolant-side heat-transfer coefficient | 57,647 W/(m²·K) |

## Thermal Resistances

| Resistance | Value | Contribution |
|---|---:|---:|
| Gas convection | 0.40558 K/W | 93.80% |
| Copper conduction | 0.00382 K/W | 0.88% |
| Coolant convection | 0.02301 K/W | 5.32% |
| **Total** | **0.43241 K/W** | **100%** |

The gas-side convection resistance dominates the Model 0 thermal-resistance network.

## Heat-Transfer Results

| Quantity | Result |
|---|---:|
| Heat-transfer rate | 1.619 kW |
| Gas-side wall heat flux | 2.58 MW/m² |
| Coolant-side wall heat flux | 2.15 MW/m² |
| Gas-side wall temperature | 343.43 K |
| Coolant-side wall temperature | 337.24 K |
| Approximate gas outlet temperature | 989.36 K |
| Approximate coolant outlet temperature | 300.78 K |

The heat fluxes differ because they are evaluated using the different inner and outer surface areas of the cylindrical wall. Both correspond to the same total heat-transfer rate.

## Coolant Boiling Screen

| Quantity | Result |
|---|---:|
| Coolant pressure used for screening | 15 bar |
| Water saturation temperature | ≈471 K |
| Calculated coolant-side wall temperature | 337.24 K |
| Saturation margin | ≈134 K |
| Model 0 outcome | No boiling predicted |

This is a saturation-temperature comparison only. It is not an onset-of-nucleate-boiling or critical-heat-flux analysis.

## Interpretation

Model 0 indicates that:

- Gas-side convection controls the calculated thermal resistance.
- Copper conduction and coolant-side convection contribute comparatively little.
- The calculated coolant-side wall temperature remains below the saturation temperature at 15 bar.
- The gas experiences only a modest bulk-temperature decrease over the 20 mm heated length.

These conclusions apply only within the documented Model 0 assumptions and correlation limitations.