#include <stdio.h>
#include <iostream>
#include <math.h>
#include "units.h"
using namespace std;

double convert_T(double T, T_unit initial, T_unit final)
{
    if (initial == final)
    {
        return T;
    }

    double c;

    switch (initial)
    {
    case C:
        return final == F ? T * 1.8 + 32 : T + 273.15;
    case F:
        c = (T - 32) / 1.8;
        return final == C ? c : c + 273.15;
    case K:
        c = T - 273.15;
        return final == C ? c : (c - 32) / 1.8;
    default:
        throw "Unrecognized unit of temperature";
    }
}

double convert_P(double P, P_unit initial, P_unit final)
{
    if (initial == final)
    {
        return P;
    }

    switch (initial)
    {
    case Pa:
        return final == bar ? P / pow(10, 5) : P / 6895;
    case bar:
        return final == Pa ? P * pow(10, 5) : P * 14.5038;
    case psi:
        return final == Pa ? P * 6895 : P / 14.5038;
    default:
        throw "Unrecognized unit of pressure";
    }
}

int main()
{
    cout << convert_P(101325, P_unit::Pa, P_unit::bar) << '\n';
    cout << convert_P(101325, P_unit::Pa, P_unit::psi) << '\n';
    cout << convert_P(1000, P_unit::psi, P_unit::bar) << '\n';
    cout << convert_P(1000, P_unit::bar, P_unit::psi) << '\n';

    return 0;
}
