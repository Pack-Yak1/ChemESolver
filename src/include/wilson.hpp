#include "units.hpp"

#ifndef wilson
#define wilson

const double IDEAL_GAS_CONST = 8.314;

class BinaryWilsonModel
{
public:
    double v1;
    double v2;
    double delta_l12;
    double delta_l21;
    double L12;
    double L21;
    bool hasTemp;

    BinaryWilsonModel();

    // Initializes an instance of BinaryWilsonModel which does not have a
    // specified temperature
    BinaryWilsonModel(double v1, double v2, double delta_l12, double delta_l21);

    // Initializes an instance of BinaryWilsonModel with a specified temperature
    // `T` in units of `t_unit`
    BinaryWilsonModel(double v1, double v2, double delta_l12, double delta_l21, double T, T_unit t_unit);

    void setT(double T, T_unit t_unit);

    // Returns the activity coefficient of component 1 for a BinaryWilsonModel
    // which already has a specified temperature
    double gamma1(double x1);

    // Destructively modifies this BinaryWilsonModel and returns the acitivity
    // coefficient of component 1 at temperature `T` in units of `t_unit`.
    double gamma1(double x1, double T, T_unit t_unit);

    double gamma2(double x2);

    // Destructively modifies this BinaryWilsonModel and returns the acitivity
    // coefficient of component 2 at temperature `T` in units of `t_unit`.
    double gamma2(double x2, double T, T_unit t_unit);
};

#endif