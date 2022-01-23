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

#endif