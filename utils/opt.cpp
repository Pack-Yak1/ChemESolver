#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "opt.h"
#include "../Thermodynamics/yx.h"

using namespace std;

opt::opt(opt_func_t f, void *context, unsigned int d)
{
    this->f = f;
    this->d = d;
    this->context = context;
}

void opt::set_polytope(vector<vector<double>> &points)
{
    this->points = points;
}

/**
 * @brief Sorts `points` in place. After sorting, `points` = `[x1, ..., xn]`
 * such that f(x1) <= ... <= f(xn).
 */
void opt::sort_by_opt_function()
{
    auto comparator = [&](vector<double> p1, vector<double> p2) -> bool
    {
        // Check if dimensions correct
        if (p1.size() != d || p2.size() != d)
        {
            throw "`sort_by_opt_function` given `points` argument containing vector of incorrect dimensions.";
        }
        double cmp = f(p1, context) - f(p2, context);
        if (cmp == 0)
        {
            // Consistent tie break needed for Nelder-Meads
            for (int i = 0; i < d; i++)
            {
                double entry_cmp = p1[i] - p2[i];
                if (entry_cmp != 0)
                {
                    return entry_cmp < 0;
                }
            }
            throw("Degenerate points compared in `sort_by_opt_function`\n");
        }
        return cmp < 0;
    };
    sort(points.begin(), points.end(), comparator);
}

vector<double> opt::get_centroid()
{
    int n = points.size();
    vector<double> output(d, 0.);
    for (int i = 0; i < n; i++)
    {
        for (int dim = 0; dim < d; dim++)
        {
            output[dim] += points[i][dim];
        }
    }
    for (int dim = 0; dim < d; dim++)
    {
        output[dim] /= n;
    }
    return output;
}

/**
 * @brief
 *
 * @param points The set of points of the convex hull
 * @param centroid The centroid of `points`
 * @param f_1 The smallest value of `f(x)` over all `x` in `points`
 * @param f_n The second largest value of `f(x)` over all `x` in `points`
 * @param f_r The value of `f(x_r)` where `x_r` is the reflection point will be stored here
 * @return The reflection point `x_r`
 */
vector<double> opt::reflect(const vector<double> &centroid, double f_1, double f_n, double &f_r)
{
    vector<double> x_r(d, 0.);
    vector<double> last = points[points.size() - 1];
    for (int dim = 0; dim < d; dim++)
    {
        x_r[dim] = centroid[dim] * (1 + ALPHA) - ALPHA * last[dim];
    }
    f_r = f(x_r, context);
    if (f_r >= f_1 && f_r < f_n)
    {
        points[points.size() - 1] = x_r;
    }
    return x_r;
}

void opt::expand(const vector<double> &centroid, const vector<double> &x_r, double f_r)
{
    vector<double> x_e(d, 0.);
    for (int dim = 0; dim < d; dim++)
    {
        x_e[dim] = centroid[dim] * (1 - BETA) + BETA * x_r[dim];
    }
    double f_e = f(x_e, context);
    points[points.size() - 1] = f_e < f_r ? x_e : x_r;
}

vector<double> opt::outside_contract(const vector<double> &centroid, const vector<double> &x_r)
{
    vector<double> x_oc(d, 0.);
    for (int dim = 0; dim < d; dim++)
    {
        x_oc[dim] = centroid[dim] * (1 - GAMMA) + GAMMA * x_r[dim];
    }
    return x_oc;
}

void opt::inside_contract(const vector<double> &centroid, const vector<double> &x_r, double f_n1)
{
    vector<double> x_ic(d, 0.);
    for (int dim = 0; dim < d; dim++)
    {
        x_ic[dim] = centroid[dim] * (1 + GAMMA) - GAMMA * x_r[dim];
    }
    double f_ic = f(x_ic, context);
    if (f_ic < f_n1)
    {
        points[points.size() - 1] = x_ic;
    }
}

void opt::shrink(const vector<double> &x_1)
{
    int n = points.size();
    for (int i = 2; i < n; i++)
    {
        for (int dim = 0; dim < d; dim++)
        {
            points[i][d] = x_1[d] * (1 - DELTA) + DELTA * points[i][d];
        }
    }
}

// Test function
double fun(const vector<double> &x, void *ctx)
{
    return -pow((x[0] - 3), 2);
}

void opt::step()
{
    double f_1 = f(points[0], context);
    double f_n = f(points[points.size() - 2], context);
    double f_n1 = f(points[points.size() - 1], context);
    vector<double> centroid = get_centroid();

    // Step 1
    sort_by_opt_function();

    // Step 2
    double f_r;
    vector<double> x_r = reflect(centroid, f_1, f_n, f_r);

    // Step 3
    if (f_r < f_1)
    {
        expand(centroid, x_r, f_r);
    }

    // Step 4
    vector<double> x_oc;
    if (f_n <= f_r && f_r < f_n1)
    {
        x_oc = outside_contract(centroid, x_r);
        double f_oc = f(x_oc, context);
        if (f_oc > f_r)
        {
            // Skip to step 6
            shrink(points[0]);
            return;
        }
        else
        {
            points[points.size() - 1] = x_oc;
        }
    }

    // Step 5
    if (f_r >= f_n1)
    {
        inside_contract(centroid, x_r, f_n1);
    }

    // Step 6
    shrink(points[0]);
}

// Optimization function. `f_data` is an instance of ModifiedRaoultModel.
// Returns the squared error of `f_data`'s pressure and the value predicted
// given `x[0]` = liquid mole fraction of component 1
// and `x[2]` = temperature.
double p_error(const std::vector<double> &x, void *f_data)
{
    ModifiedRaoultModel *m = (ModifiedRaoultModel *)f_data;
    double psat1 = m->psat_helper(m->a1, x[2]);
    double psat2 = m->psat_helper(m->a2, x[2]);
    double gamma1 = m->b.gamma1(x[0], x[2], m->t_unit);
    double gamma2 = m->b.gamma2(1 - x[0], x[2], m->t_unit);
    double p1 = x[0] * gamma1 * psat1;
    double p2 = (1 - x[0]) * gamma2 * psat2;
    return pow(m->P - p1 - p2, 2);
}

// double x_error(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
// {
//     ModifiedRaoultModel *m = (ModifiedRaoultModel *)f_data;
//     double psat1 = m->psat_helper(m->a1, x[2]);
//     double psat2 = m->psat_helper(m->a2, x[2]);
//     double gamma1 = m->b.gamma1(x[0], x[2], m->t_unit);
//     double gamma2 = m->b.gamma2(1 - x[0], x[2], m->t_unit);
//     double p1 = x[0] * gamma1 * psat1;
//     double p2 = (1 - x[0]) * gamma2 * psat2;
//     // double x1 = x[1] * m->P / gamma1 / psat1;
//     // double x2 = (1 - x[1]) * m->P / gamma2 / psat2;
//     // cout << x[0] << ',' << x[1] << ',' << x[2] << '\n';

//     return pow(x[1] - p1 / m->P, 2) + pow((1 - x[1]) - p2 / m->P, 2);
// }

int main()
{
    opt o(fun, NULL, 1);
    vector<vector<double>> v;
    for (int i = 0; i < 10; i++)
    {
        v.emplace_back(vector<double>(1, i));
    }
    o.set_polytope(v);
    // o.sort_by_opt_function();
    for (int i = 0; i < 10; i++)
    {
        cout << v[i][0] << '\n';
    }
    return 0;
}