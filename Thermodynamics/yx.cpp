#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../utils/units.h"
#include "wilson.h"
#include "antoine.h"
#include "yx.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include <fstream>

const double REL_XTOL = 1e-10;

using namespace std;

ModifiedRaoultModel::
    ModifiedRaoultModel(AntoineModel a1, AntoineModel a2, BinaryWilsonModel b,
                        T_unit t_unit, double P, P_unit p_unit)
{
    // if (a1.t_unit != t_unit)
    // {
    //     throw "The first Antoine Model provided uses a different unit of temperature than the one specified for the Modified Raoult Model.\n";
    // }
    // if (a2.t_unit != t_unit)
    // {
    //     throw "The second Antoine Model provided uses a different unit of temperature than the one specified for the Modified Raoult Model.\n";
    // }
    // if (a1.p_unit != p_unit)
    // {
    //     throw "The first Antoine Model provided uses a different unit of pressure than the one specified for the Modified Raoult Model.\n";
    // }
    // if (a2.p_unit != p_unit)
    // {
    //     throw "The second Antoine Model provided uses a different unit of pressure than the one specified for the Modified Raoult Model.\n";
    // }
    // b.setT(T, t_unit);
    this->a1 = a1;
    this->a2 = a2;
    this->b = b;
    this->t_unit = t_unit;
    this->P = P;
    this->p_unit = p_unit;

    // for (int i = 0; i <= NUM_POINTS; i++)
    // {
    //     double x1 = i / NUM_POINTS;
    //     double x2 = 1 - x1;
    //     double psat1 = a.
    // }
}

// Helper function which returns the psat of a component modeled by a1 or a2
// in units of `p_unit`, which handles unit conversion of temperature from
// `t_unit`.
double ModifiedRaoultModel::psat_helper(AntoineModel a, double T)
{
    return a.psat(convert_T(T, t_unit, a.t_unit), p_unit);
}

// Optimization function. `f_data` is an instance of ModifiedRaoultModel.
// Returns the squared error of `f_data`'s pressure and the value predicted
// given `x[0]` = mole fraction of component 1
// and `x[1]` = temperature.
double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    ModifiedRaoultModel *m = (ModifiedRaoultModel *)f_data;
    double psat1 = m->a1.psat(convert_T(x[1], m->t_unit, m->a1.t_unit), m->p_unit);
    double psat2 = m->a2.psat(convert_T(x[1], m->t_unit, m->a2.t_unit), m->p_unit);
    return pow(
        m->P - x[0] * m->b.gamma1(x[0], x[1], m->t_unit) * psat1 -
            (1 - x[0]) * m->b.gamma2(1 - x[0], x[1], m->t_unit) * psat2,
        2);
}

// Finds a value of mole fraction of component 1 which satisfies the modified
// Raoult model for temperature `T` in K.
double ModifiedRaoultModel::find_x1(double T)
{
    // x_0 refers to mole fraction of component 1, x_1 is temperature
    nlopt::opt opt(nlopt::LN_COBYLA, 2);
    std::vector<double> lower_bounds(2);
    std::vector<double> upper_bounds(2);
    ModifiedRaoultModel *f_data = this;
    lower_bounds[0] = 0;
    lower_bounds[1] = T;
    opt.set_lower_bounds(lower_bounds);
    upper_bounds[0] = 1;
    upper_bounds[1] = T;
    opt.set_upper_bounds(upper_bounds);
    opt.set_min_objective(f, f_data);
    opt.set_xtol_rel(REL_XTOL);
    std::vector<double> output(2);
    output[0] = 0.5;
    output[1] = T;
    double min_f_val;

    try
    {
        nlopt::result result = opt.optimize(output, min_f_val);
        // std::cout << "found minimum at f(" << output[0] << "," << output[1] << ") = "
        //           << std::setprecision(10) << min_f_val << std::endl;
    }
    catch (std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    return output[0];
}

// Finds a value of mole fraction of component 1 which satisfies the modified
// Raoult model for temperature `T` in units of `t_unit`.
double ModifiedRaoultModel::find_x1(double T, T_unit t_unit)
{
    return find_x1(convert_T(T, t_unit, T_unit::K));
}

double ModifiedRaoultModel::find_T(double x1)
{
    // x_0 refers to mole fraction of component 1, x_1 is temperature
    nlopt::opt opt(nlopt::LN_COBYLA, 2);
    std::vector<double> lower_bounds(2);
    std::vector<double> upper_bounds(2);
    ModifiedRaoultModel *f_data = this;
    lower_bounds[0] = x1;
    lower_bounds[1] = 0;
    opt.set_lower_bounds(lower_bounds);
    upper_bounds[0] = x1;
    upper_bounds[1] = HUGE_VAL;
    opt.set_upper_bounds(upper_bounds);
    opt.set_min_objective(f, f_data);
    opt.set_xtol_rel(REL_XTOL);
    std::vector<double> output(2);
    output[0] = x1;
    output[1] = 300;
    double min_f_val;

    try
    {
        nlopt::result result = opt.optimize(output, min_f_val);
        // std::cout << "found minimum at f(" << output[0] << "," << output[1] << ") = "
        //           << std::setprecision(10) << min_f_val << std::endl;
    }
    catch (std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    return output[1];
}

// Requires: `num_points` is no less than 2, `x_data` and `y_data` are empty
void ModifiedRaoultModel::generate_yx_data(std::vector<double> &x_data, std::vector<double> &y_data, int num_points)
{
    if (num_points < 2)
    {
        throw "`generate_yx_data` was called with `num_points` < 2.\n";
    }
    if (!x_data.empty())
    {
        throw "`generate_yx_data` was called with nonempty `x_data`.\n";
    }
    if (!y_data.empty())
    {
        throw "`generate_yx_data` was called with nonempty `y_data`.\n";
    }
    x_data.reserve(num_points);
    y_data.reserve(num_points);
    double step_size = 1. / (num_points - 1);
    ofstream output_file;
    output_file.open("out.csv");
    for (int i = 0; i < num_points; i++)
    {
        double x1 = i * step_size;
        double T = find_T(x1);
        double y1 = x1 * b.gamma1(x1, T, t_unit) * psat_helper(a1, T) / P;
        output_file << x1 << ',' << T << ',' << y1 << '\n';
        x_data.emplace_back(x1);
        y_data.emplace_back(y1);
    }
}

int main()
{
    AntoineModel a1(11.9869, 3643.32, -33.434, P_unit::bar, T_unit::K, log_type::LN);
    // std::cout << a1.tsat(329283.8093 / 100000) << '\n';
    AntoineModel a2(11.9647, 3984.93, -39.734, P_unit::bar, T_unit::K, log_type::LN);
    BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398, 363.15, T_unit::K);
    // ModifiedRaoultModel d(a, b);
    // b.setT(339.261, T_unit::K);
    // std::cout << b.gamma1(0.854518283) << '\n';
    // std::cout << b.gamma2(1 - 0.854518283) << '\n';
    ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);
    // m.find_T(0.345605278);
    // double psat1 = a1.psat(convert_T(363.15, m.t_unit, a1.t_unit), m.p_unit);
    // double psat2 = a2.psat(convert_T(363.15, m.t_unit, a2.t_unit), m.p_unit);
    // double xd = m.P - 0.060362514 * m.b.gamma1(0.060362514, 363.15, m.t_unit) * psat1 -
    //             (1 - 0.060362514) * m.b.gamma2(1 - 0.060362514, 363.15, m.t_unit) * psat2;
    std::vector<double> x_data;
    std::vector<double> y_data;

    m.generate_yx_data(x_data, y_data, 101);
    // for (double d : x_data)
    // {
    //     output_file << d << ',';
    // }
    // output_file << '\n';
    // for (double d : y_data)
    // {
    //     output_file << d << ',';
    // }
    return 0;
}