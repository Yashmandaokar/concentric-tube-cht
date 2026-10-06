# Model 0 Assumptions and Limitations

## Purpose

Model 0 is a steady, lumped, constant-property screening model. It establishes a transparent baseline for subsequent development but is not a validated design prediction.

## Assumptions

- Steady-state operation
- One-dimensional radial heat transfer
- Constant thermophysical properties
- Uniform bulk temperatures within each fluid stream
- Inlet bulk temperatures used as the thermal driving temperatures
- Smooth internal surfaces
- Single-phase liquid coolant
- Constant copper thermal conductivity
- Negligible axial conduction
- Negligible thermal radiation
- Negligible external heat loss
- Fully developed turbulent convection correlations
- Ideal-gas treatment of the gas mixture
- Ideal isentropic choking used only to estimate gas mass flow
- Uniform convective heat-transfer coefficients over the heated length

## Correlations and Applicability

The gas-side and coolant-side heat-transfer coefficients are calculated using the Gnielinski correlation with a smooth-pipe friction-factor approximation.

The coolant calculation represents turbulent annular flow using the annular hydraulic diameter.

The gas passage has a short length-to-diameter ratio:

$$
\frac{L}{D} = 2
$$

The reference gas Mach number is approximately:

$$
M_{\mathrm{ref}} \approx 0.59
$$

Consequently, the assumptions of fully developed flow and negligible compressibility are weak on the gas side. The resulting gas-side heat-transfer coefficient should be treated as an order-of-magnitude baseline.

## Gas-Flow Limitation

The ideal choking calculation determines whether choking is possible under the imposed pressure ratio and supplies a provisional gas mass flow.

It does not solve the coupled effects of:

- Wall friction
- Heat removal
- Axial pressure variation
- Axial temperature variation
- Property variation
- Compressible-flow development
- Inlet and exit losses
- Sonic-state location

Therefore, the Model 0 gas mass flow is not a prediction for the cooled, frictional, constant-area tube.

## Thermal-Model Limitation

The heat-transfer rate is based on the inlet bulk-temperature difference:

$$
\dot{Q}
=
\frac{T_{g,\mathrm{in}}-T_{c,\mathrm{in}}}
{R_g+R_{\mathrm{cond}}+R_c}
$$

The model does not account for the continuous axial evolution of fluid temperatures. The reported outlet temperatures are calculated afterward from energy balances and are approximate.

## Boiling-Screening Limitation

The coolant-side wall temperature is compared with the water saturation temperature at 15 bar.

This is a preliminary saturation-margin check only. It does not model:

- Local pressure variation
- Onset of nucleate boiling
- Subcooled boiling
- Critical heat flux
- Two-phase flow
- Flow instability

A prediction of no boiling in Model 0 therefore means only that the calculated wall temperature remains below the saturation temperature under the stated assumptions.

## Deliberately Deferred Work

The following are outside the frozen Model 0 scope:

- Entrance-effect corrections
- Correlation-sensitivity analysis
- Effectiveness–NTU formulation
- Coupled compressible momentum and energy solution
- Axially varying thermophysical properties
- Detailed boiling and critical-heat-flux correlations
- CFD validation

These effects may be introduced and documented in later model levels.