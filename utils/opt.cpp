#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "opt.h"
#include "../Thermodynamics/yx.h"
#include <float.h>

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
            throw DIMENSION_ERROR; // `sort_by_opt_function` given `points` argument containing vector of incorrect dimensions.
        }
        double f1 = f(p1, context);
        double f2 = f(p2, context);
        double cmp;
        if (isnan(f1))
        {
            cmp = 1;
        }
        else if (isnan(f2))
        {
            cmp = -1;
        }
        else
        {
            cmp = f(p1, context) - f(p2, context);
        }
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
            // throw DEGENERATE_POINTS; // "Degenerate points compared in `sort_by_opt_function`\n")
        }
        return cmp < 0;
    };
    sort(points.begin(), points.end(), comparator);
}

/**
 * @brief Find the centroid of all but the last point in `points`.
 *
 * @return vector<double>
 */
vector<double> opt::get_centroid()
{
    int n = points.size();
    vector<double> output(d, 0.);
    for (int i = 0; i < n - 1; i++)
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
 * @param x_r The value of `x_r` where `x_r` is the reflection point will be stored here
 * @return true iff the worst point in `points` was replaced with x_r.
 */
bool opt::reflect(const vector<double> &centroid, double f_1, double f_n, double &f_r, vector<double> &x_r)
{
    vector<double> last = points[points.size() - 1];
    for (int dim = 0; dim < d; dim++)
    {
        x_r[dim] = centroid[dim] * (1 + ALPHA) - ALPHA * last[dim];
    }
    f_r = f(x_r, context);
    if (f_r >= f_1 && f_r < f_n)
    {
        points[points.size() - 1] = x_r;
        return true;
    }
    return false;
}

void opt::expand(const vector<double> &centroid, const vector<double> &x_r, double f_r)
{
    vector<double> x_e(d, 0.);
    for (int dim = 0; dim < d; dim++)
    {
        x_e[dim] = centroid[dim] * (1 - BETA) + BETA * x_r[dim];
    }
    double f_e = f(x_e, context);
    cout << "Expansion point" << x_e[0] << "\n";
    if (f_e < f_r)
    {
        points[points.size() - 1] = x_e;
    }
    else
    {
        points[points.size() - 1] = x_r;
    }
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

/**
 * @brief
 *
 * @param centroid
 * @param f_n1
 * @return true iff the worst point was reflected with `x_ic` and should return to step 1
 */
bool opt::inside_contract(const vector<double> &centroid, double f_n1)
{
    vector<double> x_ic(d, 0.);
    for (int dim = 0; dim < d; dim++)
    {
        x_ic[dim] = centroid[dim] * (1 - GAMMA) + GAMMA * points[points.size() - 1][dim];
    }
    double f_ic = f(x_ic, context);
    if (f_ic < f_n1)
    {
        points[points.size() - 1] = x_ic;
        return true;
    }
    return false;
}

void opt::shrink(const vector<double> &x_1)
{
    int n = points.size();
    for (int i = 1; i < n; i++)
    {
        for (int dim = 0; dim < d; dim++)
        {
            points[i][d] = x_1[d] * (1 - DELTA) + DELTA * points[i][d];
        }
        cout << "Point " << i << " was updated to " << points[i][0] << "\n";
    }
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
    vector<double> x_r(d, 0.);
    if (reflect(centroid, f_1, f_n, f_r, x_r))
    {
        return;
    }
    f_r = isnan(f_r) ? DBL_MAX : f_r;
    cout << "Reflection point " << x_r[0] << "\n";

    // Step 3
    if (f_r < f_1)
    {
        expand(centroid, x_r, f_r);
        return;
    }

    // Step 4
    vector<double> x_oc;
    if (f_n <= f_r && f_r < f_n1)
    {
        x_oc = outside_contract(centroid, x_r);
        cout << "outside contraction point " << x_oc[0] << "\n";
        double f_oc = f(x_oc, context);
        f_oc = isnan(f_oc) ? DBL_MAX : f_oc;
        if (f_oc < f_r)
        {
            points[points.size() - 1] = x_oc;
            return;
        }
        else
        {
            // Skip to step 6
            shrink(points[0]);
            return;
        }
    }
    cout << "f(x_r) = " << f_r << ", f(x_n+1) = " << f_n1 << "\n";

    // Step 5
    if (f_r >= f_n1)
    {
        if (inside_contract(centroid, f_n1))
        {
            return;
        }
    }

    // Step 6
    shrink(points[0]);
}

void opt::set_context(void *context)
{
    this->context = context;
}

void opt::print_points()
{
    cout << "[\n";
    for (auto it = points.begin(); it != points.end(); it++)
    {
        cout << "\t[";
        for (auto v = (*it).begin(); v != (*it).end(); v++)
        {
            cout << *v << ", ";
        }
        cout << "]\n";
    }
    cout << "]\n";
}

// Test function
double fun(const vector<double> &x, void *ctx)
{
    return pow((x[0] - 3), 2);
}

typedef struct fixed_temp
{
    ModifiedRaoultModel m;
    double temp;
} * fixed_temp_t;

// Optimization function. `f_data` is an instance of ModifiedRaoultModel.
// Returns the squared error of `f_data`'s pressure and the value predicted
// given `x[0]` = liquid mole fraction of component 1
// and `x[2]` = temperature.
double p_error(const vector<double> &x, void *f_data)
{
    fixed_temp_t context = (fixed_temp_t)f_data;
    ModifiedRaoultModel *m = &(context->m);
    double temp = context->temp;
    double psat1 = m->psat_helper(m->a1, temp);
    double psat2 = m->psat_helper(m->a2, temp);
    double gamma1 = m->b.gamma1(x[0], temp, m->t_unit);
    double gamma2 = m->b.gamma2(1 - x[0], temp, m->t_unit);
    double p1 = x[0] * gamma1 * psat1;
    double p2 = (1 - x[0]) * gamma2 * psat2;
    return pow(m->P - p1 - p2, 2);
}


int main()
{
    AntoineModel a1 = ANTOINE_METHANOL;
    AntoineModel a2 = ANTOINE_WATER;
    BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398, 363.15, T_unit::K);
    ModifiedRaoultModel m(a1, a2, b, T_unit::K, 1, P_unit::bar);
    double meth_bp = a1.tsat(1);
    double water_bp = a2.tsat(1);

    double query_temperature = 350;

    // Solve for x given T
    fixed_temp_t context = (fixed_temp_t)malloc(sizeof(fixed_temp));
    context->m = m;
    context->temp = query_temperature;

    opt o(p_error, context, 1);
    vector<vector<double>> initial_points;

    for (double i = 0; i < 100; i++)
    {
        initial_points.emplace_back(vector<double>{i / 100});
    }
    // initial_points.emplace_back(vector<double>{0.5});
    // initial_points.emplace_back(vector<double>{0});
    // initial_points.emplace_back(vector<double>{1});
    // initial_points.emplace_back(vector<double>{.75});
    // initial_points.emplace_back(vector<double>{.25});
    // initial_points.emplace_back(vector<double>{1, 0, 350});
    o.set_polytope(initial_points);

    for (int i = 0; i < 1000; i++)
    {
        o.step();
        o.print_points();
    }
}

// int main()
// {

//     opt o(fun, NULL, 1);
//     vector<vector<double>> v;

//     v.emplace_back(vector<double>{0});
//     v.emplace_back(vector<double>{10});
//     v.emplace_back(vector<double>{1});
//     // for (int i = 0; i < 10; i++)
//     // {
//     //     v.emplace_back(vector<double>(1, i));
//     // }
//     o.set_polytope(v);
//     // o.sort_by_opt_function();
//     for (int i = 0; i < 10; i++)
//     {
//         o.step();
//         o.print_points();
//     }
//     return 0;
// }