#include "units.hpp"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <iostream>
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
        throw invalid_argument("Unrecognized unit of temperature.");
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
        switch (final)
        {
        case bar:
            return P / pow(10, 5);
        case psi:
            return P / 6895;
        case mmHg:
            return P / 133;
        }
        break;
    case bar:
        switch (final)
        {
        case Pa:
            return P * pow(10, 5);
        case psi:
            return P * 14.5038;
        case mmHg:
            return P * 750;
        }
        break;
    case psi:
        switch (final)
        {
        case Pa:
            return P * 6895;
        case bar:
            return P / 14.5038;
        case mmHg:
            return P * 51.715;
        }
        break;
    case mmHg:
        switch (final)
        {
        case Pa:
            return P * 133;
        case bar:
            return P / 750;
        case psi:
            return P / 51.715;
        }
        break;
    default:
        throw invalid_argument("Unrecognized unit of pressure.");
    }
    assert(0);
}

// int main()
// {
//     cout << convert_P(101325, P_unit::Pa, P_unit::bar) << '\n';
//     cout << convert_P(101325, P_unit::Pa, P_unit::psi) << '\n';
//     cout << convert_P(1000, P_unit::psi, P_unit::bar) << '\n';
//     cout << convert_P(1000, P_unit::bar, P_unit::psi) << '\n';

//     return 0;
// }
