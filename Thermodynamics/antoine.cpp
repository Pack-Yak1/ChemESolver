#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../utils/units.h"
using namespace std;

enum log_type
{
    LN,
    LOG
};

class AntoineModel
{
public:
    double A;
    double B;
    double C;
    P_unit p_unit;
    T_unit t_unit;
    log_type log_t;

    AntoineModel(double a, double b, double c, P_unit p = P_unit::Pa, T_unit t = T_unit::K, log_type l = log_type::LOG)
    {
        A = a;
        B = b;
        C = c;
        p_unit = p;
        t_unit = t;
        log_t = l;
    }

    // Returns the saturation pressure in units of `p_unit` at temeprature `T`.
    double psat(double T)
    {
        double rhs = A - B / (C + T);
        return log_t == log_type::LOG ? pow(10, rhs) : exp(rhs);
    }

    // Returns the saturation pressure in units of `p` at temeprature `T`.
    double psat(double T, P_unit p)
    {
        return convert_P(psat(T), p_unit, p);
    }

    // Returns the temperature at which the model has a saturation pressure `P`
    // in units of t_unit.
    double tsat(double P)
    {
        double lhs = log_t == log_type::LOG ? log10(P) : log(P);
        return B / (A - lhs) - C;
    }

    // Returns the temperature at which the model has a saturation pressure `P`
    // in units of t.
    double tsat(double P, T_unit t)
    {
        return convert_T(tsat(P), t_unit, t);
    }
};

int main()
{
    AntoineModel a(11.9869, 3643.32, -33.434, P_unit::bar, T_unit::K, log_type::LN);
    cout << a.tsat(329283.8093 / 100000) << '\n';
    return 0;
}
