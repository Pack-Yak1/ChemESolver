# ChemESolver

Calculation package for Chemical Engineering problems

## Distillation

McCabe-Thiele Analysis for non-ideal liquid mixtures. Currently supports the
modified Raoult's Law model for vapor-liquid equilibrium. Support for ideal
Raoult's Law in the pipeline. 

## Thermodynamics 

Models for determining physical properties of substances/mixtures. E.g. boiling
point, activity coefficients, vapor pressure. Support for additional models as 
well as a unified interface in the pipeline.

### Features:

#### Distillation

1. Determination of minimum reflux ratio under optionally specified feed conditions
2. Calculation of operating lines given flow rates or ratios and desired product purities
3. Calculation of pseudoequilibrium curve under specified operating conditions and Murphree efficiency
4. Determination of stage counts under specified conditions
5. Calculation of tray compositions for McCabe-Thiele diagram plotting

#### Thermodynamics

1. Antoine Model for calculation of vapor pressures/boiling points
2. Binary Wilson Model for calculation of liquid phase activity coefficients
3. VLE models utilizing Antoine and Binary Wilson Models to determine binary mixture properties and solve for equilibrium temperature/compositions

#### Misc. (utils)

1. Support for implicit unit conversions
2. Rudimentary vector operations (E.g. arithmetic, display, linspace/mesh generation)
3. Solver using Nelder-Meads algorithm for non-gradient based optimization