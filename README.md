# ChemESolver

Calculation package for Chemical Engineering problems

## Usage

Header files specifying interfaces are in `src/include`. Source files for 
ChemESolver's modules are in `src/lib`. 

To include ChemESolver's modules in your application, simply run `make all`. 
This produces shared and static libraries in `bin/shared` and `bin/static`, 
respectively. To use the static library, compile your program with the flags

```-L<path>bin/static -I<path>src/include -lChemESolver```

To use the shared library, compile your program with the flags

```-L<path>bin/shared -I<path>src/include -lChemESolver -Wl,-rpath=<path>bin/shared```

after your program's source file(s).

Alternatively, users can write C++ programs in `src/apps` and run `make all` to
compile their programs with all dependencies from this library linked. 
Executables for each C++ source file will be stored in `bin`.

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
4. Equilibrium diagram (Txy, xy) diagram data computation

#### Misc. (utils)

1. Support for implicit unit conversions
2. Rudimentary vector operations (E.g. arithmetic, display, linspace/mesh generation)
3. Solver using Adaptive Nelder-Mead Simplex algorithm for non-gradient based optimization