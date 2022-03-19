#include "../utils/units.h"

#ifndef antoine
#define antoine

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

    AntoineModel();

    AntoineModel(double a, double b, double c, P_unit p = P_unit::Pa, T_unit t = T_unit::K, log_type l = log_type::LOG);

    // Returns the saturation pressure in units of `p_unit` at temeprature `T`.
    double psat(double T);

    // Returns the saturation pressure in units of `p` at temeprature `T`.
    double psat(double T, P_unit p);

    // Returns the temperature at which the model has a saturation pressure `P`
    // in units of t_unit.
    double tsat(double P);

    // Returns the temperature at which the model has a saturation pressure `P`
    // in units of t.
    double tsat(double P, T_unit t);
};

// Antoine model for water with units of `bar`, `K`.
const AntoineModel ANTOINE_WATER =
    AntoineModel(
        11.9647, 3984.93, -39.734, P_unit::bar, T_unit::K, log_type::LN);

// Antoine model for methanol with units of `bar`, `K`.
const AntoineModel ANTOINE_METHANOL =
    AntoineModel(
        11.9869, 3643.32, -33.434, P_unit::bar, T_unit::K, log_type::LN);

#endif