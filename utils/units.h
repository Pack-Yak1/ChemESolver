#ifndef units
#define units

enum T_unit
{
    C,
    F,
    K
};

enum P_unit
{
    Pa,
    bar,
    psi
};

double convert_T(double T, T_unit initial, T_unit final);

double convert_P(double P, P_unit initial, P_unit final);

#endif