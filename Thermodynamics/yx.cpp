#include "yx.h"

#include <math.h>
#include <stdio.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>

#include "../utils/coords.h"
#include "../utils/opt.h"
#include "../utils/units.h"
#include "antoine.h"
#include "gnuplot-iostream.h"
#include "wilson.h"

using namespace std;

#define x_var 0
#define y_var 1
#define T_var 2

class opt_spec
{
public:
    ModifiedRaoultModel *m;
    double val;
    int constraint;

    opt_spec(ModifiedRaoultModel *m, double val, int constraint)
    {
        this->m = m;
        this->val = val;
        // 0 for x1, 1 for y1, 2 for T.
        this->constraint = constraint;
    }
};

ModifiedRaoultModel::ModifiedRaoultModel()
{
    this->a1 = AntoineModel();
    this->a2 = AntoineModel();
    this->b = BinaryWilsonModel();
    this->t_unit = T_unit::K;
    this->P = 101325;
    this->p_unit = P_unit::Pa;
}

ModifiedRaoultModel::ModifiedRaoultModel(const ModifiedRaoultModel &m)
{
    this->a1 = m.a1;
    this->a2 = m.a2;
    this->b = m.b;
    this->t_unit = m.t_unit;
    this->P = m.P;
    this->p_unit = m.p_unit;
}

ModifiedRaoultModel::ModifiedRaoultModel(AntoineModel a1, AntoineModel a2,
                                         BinaryWilsonModel b, T_unit t_unit,
                                         double P, P_unit p_unit)
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

/**
 * @brief Loss function. Squared error of predicted vs actual pressure of the
 * system. If context constrains x1 or T, then `x` is the other quantity.
 *
 * @param x `x[0]` = liquid phase mole fraction of component 1. Must be size 1.
 * @param context points to an opt_spec object. `val` must be the fixed value of
 * temperature in the system
 * @return The value of the loss function.
 */
double p_error(const vector<double> &x, void *context)
{
    assert(x.size() == 1);
    opt_spec *spec = (opt_spec *)context;
    ModifiedRaoultModel *m = spec->m;
    int constraint = spec->constraint;
    double liq;
    double temp;

    switch (constraint)
    {
    case x_var:
        liq = spec->val;
        temp = x[0];
        break;
    case T_var:
        temp = spec->val;
        liq = x[0];
        break;
    default:
        assert(0);
    }
    double psat1 = m->psat_helper(m->a1, temp);
    double psat2 = m->psat_helper(m->a2, temp);
    double gamma1 = m->b.gamma1(liq, temp, m->t_unit);
    double gamma2 = m->b.gamma2(1 - liq, temp, m->t_unit);
    double p1 = liq * gamma1 * psat1;
    double p2 = (1 - liq) * gamma2 * psat2;
    return pow(m->P - p1 - p2, 2);
}

/**
 * @brief
 *
 * @param x
 * @param context
 * @return double
 */
double x_error(const std::vector<double> &x, void *context)
{
    assert(x.size() == 2);
    opt_spec *spec = (opt_spec *)context;
    ModifiedRaoultModel *m = spec->m;
    int constraint = spec->constraint;
    switch (constraint)
    {
    case y_var:
        break;
    default:
        assert(0);
    }
    double liq = x[0];
    double temp = x[1];
    double vap = spec->val;

    double psat1 = m->psat_helper(m->a1, temp);
    double psat2 = m->psat_helper(m->a2, temp);
    double gamma1 = m->b.gamma1(liq, temp, m->t_unit);
    double gamma2 = m->b.gamma2(1 - liq, temp, m->t_unit);
    double p1 = liq * gamma1 * psat1;
    double p2 = (1 - liq) * gamma2 * psat2;

    return pow(vap - p1 / m->P, 2) + pow((1 - vap) - p2 / m->P, 2);
}

// Finds a value of mole fraction of component 1 which satisfies the modified
// Raoult model for temperature `T` in K.
vector<double> ModifiedRaoultModel::solve_from_T(double T)
{
    opt_spec context(this, T, T_var); // Specify that we fix the value of T
    opt solver(p_error, &context, 1);
    vector<double> lb{0.}; // Solution is bounded by [0, 1]
    vector<double> ub{1.};

    solution *soln = solver.solve(lb, ub);
    double x1_soln = soln->x[0];
    double psat1 = psat_helper(a1, T);
    double gamma1 = b.gamma1(x1_soln, T, t_unit);
    double p1 = x1_soln * gamma1 * psat1;
    delete soln;
    return vector<double>{x1_soln, p1 / P, T};
}

// Finds a value of mole fraction of component 1 which satisfies the modified
// Raoult model for temperature `T` in units of `t_unit`.
vector<double> ModifiedRaoultModel::solve_from_T(double T, T_unit t_unit)
{
    return solve_from_T(convert_T(T, t_unit, T_unit::K));
}

vector<double> ModifiedRaoultModel::solve_from_y1(double y1)
{
    opt_spec context(this, y1, y_var); // Specify that we fix the value of x1
    opt solver(x_error, &context, 2);
    double bp1 = tsat_helper(a1, P);
    double bp2 = tsat_helper(a2, P);
    // Solution x1 bounded by [0, 1], T bounded by bp's.
    vector<double> lb{0, 0};
    vector<double> ub{1, 400};

    solution *soln = solver.solve(lb, ub);
    double x = soln->x[0];
    double t = soln->x[1];
    delete soln;
    return vector<double>{x, y1, t};
}

vector<double> ModifiedRaoultModel::solve_from_x1(double x1)
{
    opt_spec context(this, x1, x_var); // Specify that we fix the value of x1
    opt solver(p_error, &context, 1);
    double bp1 = tsat_helper(a1, P);
    double bp2 = tsat_helper(a2, P);
    vector<double> lb{min(bp1, bp2)}; // Lower boiling point must be the lb
    vector<double> ub{max(bp1, bp2)}; // Higher boiling point must be the ub

    solution *soln = solver.solve(lb, ub);
    double t_soln = soln->x[0];
    double psat1 = psat_helper(a1, t_soln);
    double gamma1 = b.gamma1(x1, t_soln, t_unit);
    double p1 = x1 * gamma1 * psat1;
    delete soln;
    return vector<double>{x1, p1 / P, t_soln};
}

void ModifiedRaoultModel::set_Ty(double x1, double &y1, double &T)
{
    vector<double> answer = solve_from_x1(x1);
    y1 = answer[1];
    T = answer[2];
}

void ModifiedRaoultModel::set_Tx(double &x1, double y1, double &T)
{
    vector<double> answer = solve_from_y1(y1);
    x1 = answer[0];
    T = answer[2];
}

void ModifiedRaoultModel::set_xy(double &x1, double &y1, double T)
{
    vector<double> answer = solve_from_T(T);
    x1 = answer[0];
    y1 = answer[1];
}

// Throws exceptions parametrized on `fn_name` if num_points is less than 2,
// or if start and end are not within [0, 1], or if start is greater than end.
void check_bounds(int num_points, double start, double end, string fn_name)
{
    if (num_points < 2)
    {
        throw invalid_argument("`" + fn_name +
                               "` was called with `num_points` < 2.");
    }
    if (start < 0 || end > 1)
    {
        throw invalid_argument("`" + fn_name +
                               "` was called with `start` < 0 or `end` > 1.");
    }
    if (start > end)
    {
        throw invalid_argument("`" + fn_name +
                               "` was called with `start` > `end`.");
    }
}

void ModifiedRaoultModel::generate_Txy_single_thread(double start,
                                                     double step_size,
                                                     double num_steps, int pos,
                                                     vector<Point> &Tx_data,
                                                     vector<Point> &Ty_data,
                                                     mutex &mtx)
{
    double y1, T;
    vector<Point> thread_Tx;
    vector<Point> thread_Ty;
    for (int i = 0; i < num_steps; i++)
    {
        double x1 = start + i * step_size;
        set_Ty(x1, y1, T);
        thread_Tx.emplace_back(x1, T);
        thread_Ty.emplace_back(y1, T);
    }

    // mtx.lock();
    for (int i = 0; i < num_steps; i++)
    {
        Tx_data[pos + i] = thread_Tx[i];
        Ty_data[pos + i] = thread_Ty[i];
    }
    // mtx.unlock();
}

void ModifiedRaoultModel::generate_Txy_data(int num_points,
                                            vector<Point> &Tx_data,
                                            vector<Point> &Ty_data,
                                            double start, double end,
                                            int n_workers)
{
    check_bounds(num_points, start, end, "generate_Txy_data");
    if (!Tx_data.empty())
    {
        throw invalid_argument(
            "`generate_yx_data` was called with nonempty `Tx_data`.");
    }
    if (!Ty_data.empty())
    {
        throw invalid_argument(
            "`generate_yx_data` was called with nonempty `Ty_data`.");
    }
    Tx_data.reserve(num_points);
    Ty_data.reserve(num_points);
    double step_size = (end - start) / (num_points - 1);
    double x1, y1, T;

    // Mutual exclusion actually isn't needed here. We just don't want cores to
    // keep prefetching the entire output vector to their caches and invalidate
    // other cores' cache lines.
    mutex answer_mtx;

    // Assign each worker a partition of output graph
    vector<thread> workers;
    for (int i = 0; i < n_workers; i++)
    {
        int pos = num_points * i / n_workers;
        int num_steps = num_points * (i + 1) / n_workers - pos;
        double thread_start = pos * step_size + start;
        workers.emplace_back(
            thread(&ModifiedRaoultModel::generate_Txy_single_thread,
                   ModifiedRaoultModel(*this), thread_start, step_size,
                   num_steps, pos, std::ref(Tx_data), std::ref(Ty_data),
                   std::ref(answer_mtx)));
    }

    // Await completion of all workers
    for (auto &th : workers)
    {
        th.join();
    }
}

void ModifiedRaoultModel::write_Txy_data(int num_points,
                                         ostream &o,
                                         string delim,
                                         string line_break,
                                         double start,
                                         double end,
                                         int n_workers)
{
    check_bounds(num_points, start, end, "write_Txy_data");
    vector<Point> Tx_data;
    vector<Point> Ty_data;
    generate_Txy_data(num_points, Tx_data, Ty_data, start, end, n_workers);
    for (int i = 0; i < num_points; i++)
    {
        o << Tx_data[i].x << delim << Ty_data[i].x << delim << Tx_data[i].y << line_break;
    }
}
