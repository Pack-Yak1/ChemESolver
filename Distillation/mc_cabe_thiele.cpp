#include "../Thermodynamics/antoine.h"
#include "../Thermodynamics/wilson.h"
#include "../Thermodynamics/yx.h"
#include "../utils/units.h"
#include "../utils/point.h"
#include "mc_cabe_thiele.h"
#include <nlopt.hpp>
#include <vector>

using namespace std;

MT::MT(ModifiedRaoultModel m, double xD, double xB, double xF, bool total_reflux)
{
    this->m = m;
    this->xD = xD;
    this->xB = xB;
    this->xF = xF;
    this->total_reflux = total_reflux;
}

vector<double> MT::rectifying_line(double R)
{
    vector<double> output;
    output.reserve(2);
    output.emplace_back(R / (R + 1));
    output.emplace_back(xD / (R + 1));
    // cout << output[0] << "," << output[1] << "\n";
    return output;
}

vector<double> MT::rectifying_line(double L, double D)
{
    double R = L / D;
    return rectifying_line(R);
}

vector<double> MT::stripping_line(double VB)
{
    vector<double> output;
    output.reserve(2);
    output.emplace_back((VB + 1) / VB);
    output.emplace_back(-xB / VB);
    // cout << output[0] << "," << output[1] << "\n";
    return output;
}

vector<double> MT::stripping_line(double Vbar, double B)
{
    double VB = Vbar / B;
    return stripping_line(VB);
}

// Assumes saturated liquid feed
double MT::min_reflux()
{
    double yF = xF;
    double T = 300;
    m.find_Ty(xF, yF, T);
    // cout << xF << ',' << yF << ',' << T << '\n';
    double gradient = Line(xF, yF, xD, xD).gradient();
    return -gradient / (gradient - 1);
}

// x_0 is liq mole fraction, x_1 is T, x_2 is vap mole fraction
double q_line_constraint(const std::vector<double> &x, std::vector<double> &grad, void *constraint_data)
{
    vector<double> *v = (vector<double> *)constraint_data;
    double q_line_grad = (*v)[0];
    double q_line_intersect = (*v)[1];
    // cout << x[0] << ',' << x[1] << ',' << x[2] << '\n';
    // cout << q_line_grad << ',' << q_line_intersect << '\n';
    // cout << "constraint error: " << x[1] - q_line_grad * x[0] - q_line_intersect << "\n\n";
    return x[1] - q_line_grad * x[0] - q_line_intersect;
}

double y_constraint(const std::vector<double> &x, std::vector<double> &grad, void *constraint_data)
{
    ModifiedRaoultModel *m = (ModifiedRaoultModel *)constraint_data;
    double gamma = m->b.gamma1(x[0], x[2], m->t_unit);
    double psat = m->psat_helper(m->a1, x[2]);
    // double psat = m->a1.psat(convert_T(x[2], m->t_unit, m->a1.t_unit), m->p_unit);
    return gamma * x[0] * psat / m->P - x[1];
}

// Finds minimum reflux ratio for a specified q-value
double MT::min_reflux(double q)
{
    if (q == 1)
    {
        return min_reflux();
    }
    vector<double> constraint_data;
    constraint_data.reserve(2);
    constraint_data.emplace_back(-q / (1. - q));
    constraint_data.emplace_back(-xF / (q - 1.));

    // x_0, x_1 is liq, vap mole fraction, x_2 is T
    nlopt::opt opt(nlopt::LN_COBYLA, 3);
    std::vector<double> lower_bounds(3);
    std::vector<double> upper_bounds(3);
    ModifiedRaoultModel *f_data = &m;
    double boiling_point_1 = m.tsat_helper(m.a1, m.P);
    double boiling_point_2 = m.tsat_helper(m.a2, m.P);
    lower_bounds[0] = 0;
    lower_bounds[1] = 0;
    lower_bounds[2] = min(boiling_point_1, boiling_point_2);
    opt.set_lower_bounds(lower_bounds);
    upper_bounds[0] = 1;
    upper_bounds[1] = 1;
    upper_bounds[2] = max(boiling_point_1, boiling_point_2);
    opt.set_upper_bounds(upper_bounds);
    opt.set_min_objective(yx::f, f_data);
    opt.add_equality_constraint(q_line_constraint, &constraint_data, 1e-8);
    opt.add_equality_constraint(y_constraint, &m, 1e-8);
    opt.set_xtol_rel(yx::REL_XTOL);
    std::vector<double> output(3);
    output[0] = xF;
    output[1] = 1;
    output[2] = lower_bounds[2];
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
    cout << output[0] << ',' << output[1] << ',' << output[2] << '\n';
    double gradient = Line(output[0], output[1], xD, xD).gradient();
    return -gradient / (gradient - 1);
}

vector<Point> MT::pseudo_equilibrium_curve(int num_points, double efficiency, vector<double> rect_line, vector<double> strip_line)
{
    Point intersection = Line::intersection(rect_line, strip_line);
    vector<Point> Tx_data;
    vector<Point> Ty_data;
    m.generate_Txy_data(num_points, Tx_data, Ty_data, xB, xD);
    vector<Point> output;
    output.reserve(num_points);
    double step_size = (xD - xB) / (num_points - 1);
    for (int i = 0; i < num_points; i++)
    {
        double x = Tx_data[i].x;
        double y = Ty_data[i].x;
        double y_op;
        if (total_reflux)
        {
            y_op = x;
        }
        else
        {
            if (x < intersection.x)
            {
                y_op = strip_line[0] * x + strip_line[1];
            }
            else
            {
                y_op = rect_line[0] * x + rect_line[1];
            }
        }
        double y_out = (y - y_op) * efficiency + y_op;
        output.emplace_back(x, y_out);
    }
    return output;
    // cout << intersection << "\n";
    // return p;
}

vector<Point> MT::pseudo_equilibrium_curve(int num_points, double efficiency,
                                           double R, double VB)
{
    return pseudo_equilibrium_curve(
        num_points, efficiency, rectifying_line(R), stripping_line(VB));
}

int main()
{
    // AntoineModel a1 = antoine::ANTOINE_METHANOL;
    // AntoineModel a2 = antoine::ANTOINE_WATER;
    // cout << a2.tsat(99250 / 1e5) << '\n';

    // BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398, 363.15, T_unit::K);
    // ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);
    // MT mt(m, 0.91, 0.002895894, 0.123287671, false);
    // cout << mt.min_reflux(1.121256386) << '\n';
    // vector<Point> peq = mt.pseudo_equilibrium_curve(101, 0.1, 1.86805487, 0.829084445);

    // ofstream output_file;
    // output_file.open("test_peq.csv");
    // // mt.m.write_Txy_data(101, output_file);
    // for (int i = 0; i < peq.size(); i++)
    // {
    //     output_file << peq[i].x << "," << peq[i].y << '\n';
    // }
    return 0;
}