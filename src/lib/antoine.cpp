#include <stdio.h>
#include <iostream>
#include <math.h>
#include "units.h"
#include "antoine.h"
using namespace std;

AntoineModel::AntoineModel()
{
    A = 0;
    B = 0;
    C = 0;
    // TODO: type rest of defaults
}

AntoineModel::AntoineModel(double a, double b, double c, P_unit p,
                           T_unit t, log_type l)
{
    A = a;
    B = b;
    C = c;
    p_unit = p;
    t_unit = t;
    log_t = l;
}

// Returns the saturation pressure in units of `p_unit` at temeprature `T`.
double AntoineModel::psat(double T)
{
    double rhs = A - B / (C + T);
    return log_t == log_type::LOG ? pow(10, rhs) : exp(rhs);
}

// Returns the saturation pressure in units of `p` at temeprature `T`.
double AntoineModel::psat(double T, P_unit p)
{
    return convert_P(psat(T), p_unit, p);
}

// Returns the temperature at which the model has a saturation pressure `P`
// in units of t_unit.
double AntoineModel::tsat(double P)
{
    double lhs = log_t == log_type::LOG ? log10(P) : log(P);
    return B / (A - lhs) - C;
}

// Returns the temperature at which the model has a saturation pressure `P`
// in units of t.
double AntoineModel::tsat(double P, T_unit t)
{
    return convert_T(tsat(P), t_unit, t);
}

// int main()
// {

// }
