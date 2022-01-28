#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../utils/units.h"
#include "../utils/point.h"
#include "wilson.h"
#include "antoine.h"
#include "yx.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include <fstream>
#include <cmath>
#include "gnuplot-iostream.h"
// #include <boost/tuple/tuple.hpp>
// #include "gnuplot-iostream.h"

using namespace std;

ModifiedRaoultModel::ModifiedRaoultModel()
{
    this->a1 = AntoineModel();
    this->a2 = AntoineModel();
    this->b = BinaryWilsonModel();
    this->t_unit = T_unit::K;
    this->P = 101325;
    this->p_unit = P_unit::Pa;
}

ModifiedRaoultModel::
    ModifiedRaoultModel(AntoineModel a1, AntoineModel a2, BinaryWilsonModel b,
                        T_unit t_unit, double P, P_unit p_unit)
{
    this->a1 = a1;
    this->a2 = a2;
    this->b = b;
    this->t_unit = t_unit;
    this->P = P;
    this->p_unit = p_unit;
}

double ModifiedRaoultModel::psat_helper(AntoineModel a, double T)
{
    return a.psat(convert_T(T, t_unit, a.t_unit), p_unit);
}

double ModifiedRaoultModel::tsat_helper(AntoineModel a, double P)
{
    return a.tsat(convert_P(P, p_unit, a.p_unit), t_unit);
}

// Optimization function. `f_data` is an instance of ModifiedRaoultModel.
// Returns the squared error of `f_data`'s pressure and the value predicted
// given `x[0]` = liquid mole fraction of component 1
// and `x[2]` = temperature.
double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    ModifiedRaoultModel *m = (ModifiedRaoultModel *)f_data;
    double psat1 = m->psat_helper(m->a1, x[2]);
    double psat2 = m->psat_helper(m->a2, x[2]);
    double gamma1 = m->b.gamma1(x[0], x[2], m->t_unit);
    double gamma2 = m->b.gamma2(1 - x[0], x[2], m->t_unit);
    double p1 = x[0] * gamma1 * psat1;
    double p2 = (1 - x[0]) * gamma2 * psat2;
    double T_err = pow(m->P - p1 - p2, 2);
    return T_err;
}

// Finds a value of mole fraction of component 1 which satisfies the modified
// Raoult model for temperature `T` in K.
double ModifiedRaoultModel::find_x1(double T)
{
    // x_0 refers to mole fraction of component 1, x_1 is temperature
    nlopt::opt opt(nlopt::LN_COBYLA, 3);
    std::vector<double> lower_bounds(3);
    std::vector<double> upper_bounds(3);
    ModifiedRaoultModel *f_data = this;
    lower_bounds[0] = 0;
    lower_bounds[1] = 0;
    lower_bounds[2] = T;
    opt.set_lower_bounds(lower_bounds);
    upper_bounds[0] = 1;
    upper_bounds[1] = 1;
    upper_bounds[2] = T;
    opt.set_upper_bounds(upper_bounds);
    opt.set_min_objective(f, f_data);
    opt.set_xtol_rel(REL_XTOL);
    std::vector<double> output(3);
    output[0] = 0.5;
    output[1] = 0.5;
    output[2] = T;
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
    // x_0, x_1 is liq, vap mole fraction, x_2 is T
    nlopt::opt opt(nlopt::LN_COBYLA, 3);
    std::vector<double> lower_bounds(3);
    std::vector<double> upper_bounds(3);
    ModifiedRaoultModel *f_data = this;
    lower_bounds[0] = x1;
    lower_bounds[1] = 0;
    lower_bounds[2] = 0;
    opt.set_lower_bounds(lower_bounds);
    upper_bounds[0] = x1;
    upper_bounds[1] = 1;
    upper_bounds[2] = HUGE_VAL;
    opt.set_upper_bounds(upper_bounds);
    opt.set_min_objective(f, f_data);
    opt.set_xtol_rel(REL_XTOL);
    std::vector<double> output(3);
    output[0] = x1;
    output[1] = x1;
    output[2] = 300;
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
    return output[2];
}

void ModifiedRaoultModel::find_Ty(double x1, double &y1, double &T)
{
    T = find_T(x1);
    y1 = x1 * b.gamma1(x1, T, t_unit) * psat_helper(a1, T) / P;
}

// Throws exceptions parametrized on `fn_name` if num_points is less than 2,
// or if start and end are not within [0, 1], or if start is greater than end.
void check_bounds(int num_points, double start, double end, string fn_name)
{
    if (num_points < 2)
    {
        throw "`" + fn_name + "` was called with `num_points` < 2.\n";
    }
    if (start < 0 || end > 1)
    {
        throw "`" + fn_name + "` was called with `start` < 0 or `end` > 1.\n";
    }
    if (start > end)
    {
        throw "`" + fn_name + "` was called with `start` > `end`.\n";
    }
}

void ModifiedRaoultModel::
    generate_Txy_data(int num_points, vector<Point> &Tx_data,
                      vector<Point> &Ty_data, double start, double end)
{
    check_bounds(num_points, start, end, "generate_Txy_data");
    if (!Tx_data.empty())
    {
        throw "`generate_yx_data` was called with nonempty `Tx_data`.\n";
    }
    if (!Ty_data.empty())
    {
        throw "`generate_yx_data` was called with nonempty `Ty_data`.\n";
    }
    Tx_data.reserve(num_points);
    Ty_data.reserve(num_points);
    double step_size = (end - start) / (num_points - 1);
    double x1, y1, T;
    for (int i = 0; i < num_points; i++)
    {
        x1 = i * step_size + start;
        find_Ty(x1, y1, T);
        Tx_data.emplace_back(x1, T);
        Ty_data.emplace_back(y1, T);
    }
}

void ModifiedRaoultModel::write_Txy_data(int num_points, ostream &o,
                                         string delim, string line_break,
                                         double start, double end)
{
    check_bounds(num_points, start, end, "write_Txy_data");

    double step_size = 1. / (num_points - 1);
    double x1, y1, T;
    for (int i = 0; i < num_points; i++)
    {
        x1 = i * step_size;
        find_Ty(x1, y1, T);
        o << x1 << delim << y1 << delim << T << line_break;
    }
}

// int main()
// {
//     AntoineModel a1(11.9869, 3643.32, -33.434, P_unit::bar, T_unit::K, log_type::LN);
//     // std::cout << a1.tsat(329283.8093 / 100000) << '\n';
//     AntoineModel a2(11.9647, 3984.93, -39.734, P_unit::bar, T_unit::K, log_type::LN);
//     BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398, 363.15, T_unit::K);
//     // ModifiedRaoultModel d(a, b);
//     // b.setT(339.261, T_unit::K);
//     // std::cout << b.gamma1(0.854518283) << '\n';
//     // std::cout << b.gamma2(1 - 0.854518283) << '\n';
//     ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);
//     cout << m.find_x1(363.15) << '\n';
//     // double psat1 = a1.psat(convert_T(363.15, m.t_unit, a1.t_unit), m.p_unit);
//     // double psat2 = a2.psat(convert_T(363.15, m.t_unit, a2.t_unit), m.p_unit);
//     // double xd = m.P - 0.060362514 * m.b.gamma1(0.060362514, 363.15, m.t_unit) * psat1 -
//     //             (1 - 0.060362514) * m.b.gamma2(1 - 0.060362514, 363.15, m.t_unit) * psat2;
//     // std::vector<Point> Tx_data;
//     // std::vector<Point> Ty_data;

//     // ofstream output_file;
//     // output_file.open("out.dat");

//     // m.generate_Txy_data(101, Tx_data, Ty_data);
//     // for (Point p : Tx_data)
//     // {
//     //     cout << p << '\n';
//     // }
//     // cout << '\n';
//     // for (Point p : Ty_data)
//     // {
//     //     cout << p << '\n';
//     // }
//     // cout << '\n';
//     // m.write_Txy_data(101, output_file);
//     // Gnuplot gp;
//     return 0;
// }
