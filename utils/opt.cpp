#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "opt.h"
#include "../Thermodynamics/yx.h" // TODO: remove after testing
#include <float.h>
#include "vector_ops.h"

using namespace std;

solution::solution(vector<double> x, double fx)
{
    this->fx = fx;
    this->x = x;
}

opt::opt(opt_func_t f, void *context, unsigned int d)
{
    this->f = f;
    this->d = d;
    this->context = context;
    this->last_stddev = 0;
    this->points = vector<vector<double>>();
}

void opt::set_polytope(vector<vector<double>> &points)
{
    this->points = points;
}

void opt::set_context(void *context)
{
    this->context = context;
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
            // `sort_by_opt_function` given `points` argument containing vector
            // of incorrect dimensions.
            throw invalid_argument("Dimensions of vectors do not match in `sort_by_opt_function`. Ensure that all points provided to `set_polytope` match the dimension of the problem.");
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
        }
        return cmp < 0;
    };
    sort(points.begin(), points.end(), comparator);
}

/**
 * @brief Find the centroid of all but the last point in `points`.
 *
 * @return vector<double> The centroid of all but the worst (highest f(x)) point
 * in `points`
 */
vector<double> opt::get_centroid()
{
    int n = points.size();
    vector<double> output = multiply(1 / n, sum(points));
    return output;
}

/**
 * @brief Carries out the reflection step of the Nelder-Mead algorithm in place,
 * and store the value of the reflection point in `x_r`.
 *
 * @param centroid The centroid of all but the worst point in `points`
 * @param f_1 The smallest value of `f(x)` over all `x` in `points`
 * @param f_n The second largest value of `f(x)` over all `x` in `points`
 * @param f_r The value of `f(x_r)` where `x_r` is the reflection point will be
 * stored in the address of this variable
 * @param x_r The value of `x_r` where `x_r` is the reflection point will be
 * stored in the address of this variable. Must be initialized as d-dimensional
 * zero vector
 * @return true iff the worst point in `points` was replaced with x_r.
 */
bool opt::reflect(const vector<double> &centroid, double f_1, double f_n, double &f_r, vector<double> &x_r)
{
    vector<double> last = points[points.size() - 1];
    sum_in_place(x_r, multiply(1 + ALPHA, centroid));
    sum_in_place(x_r, multiply(-1, last));
    f_r = f(x_r, context);
    if (f_r >= f_1 && f_r < f_n)
    {
        points[points.size() - 1] = x_r;
        return true;
    }
    return false;
}

/**
 * @brief Carries out the expansion step of the Nelder-Mead algorithm in place
 *
 * @param centroid The centroid of all but the worst point in `points`
 * @param x_r The reflection point of `points`, as calculated by `reflect`
 * @param f_r The value of the objective function at the reflection point.
 */
void opt::expand(const vector<double> &centroid, const vector<double> &x_r, double f_r)
{
    vector<double> x_e(d, 0.);
    sum_in_place(x_e, multiply(1 - BETA, centroid));
    sum_in_place(x_e, multiply(BETA, x_r));
    double f_e = f(x_e, context);
#ifdef DEBUG
    cout << "Expansion point ";
    vector_println(x_e);
#endif
    if (f_e < f_r)
    {
        points[points.size() - 1] = x_e;
    }
    else
    {
        points[points.size() - 1] = x_r;
    }
}

/**
 * @brief Calculates the external contraction point for the Nelder-Mead
 * algorithm.
 *
 * @param centroid The centroid of all but the worst point in `points`
 * @param x_r The reflection point of `points`, as calculated by `reflect`
 * @return vector<double> The external contraction point of `points`
 */
vector<double> opt::outside_contract(const vector<double> &centroid, const vector<double> &x_r)
{
    vector<double> x_oc(d, 0.);
    sum_in_place(x_oc, multiply(1 - GAMMA, centroid));
    sum_in_place(x_oc, multiply(GAMMA, x_r));
    return x_oc;
}

/**
 * @brief Carries out the internal contraction step of Nelder-Mead algorithm in
 * place
 *
 * @param centroid The centroid of all but the worst point in `points`
 * @param f_n1 The value of the objective function at the worst (highest) point
 * in `points`
 * @return true iff the worst point was reflected with `x_ic` and should return
 * to step 1
 */
bool opt::inside_contract(const vector<double> &centroid, double f_n1)
{
    vector<double> x_ic(d, 0.);
    sum_in_place(x_ic, multiply(1 - GAMMA, centroid));
    sum_in_place(x_ic, multiply(GAMMA, points.back()));
    double f_ic = f(x_ic, context);
    if (f_ic < f_n1)
    {
        points[points.size() - 1] = x_ic;
        return true;
    }
    return false;
}

/**
 * @brief Carries out the shrink step of the Nelder-Mead algorithm in place.
 *
 * @param x_1 The best point at the current step of the algorithm, i.e. f(x_1)
 * where x_1 is the first point returned by `sort_by_opt_function`.
 */
void opt::shrink(const vector<double> &x_1)
{
    int n = points.size();
    for (int i = 1; i < n; i++)
    {
        multiply_in_place(DELTA, points[i]);
        sum_in_place(points[i], multiply(1 - DELTA, x_1));
#ifdef DEBUG
        cout << "Point " << i << " was updated to " << points[i][0] << "\n";
#endif
    }
}

/**
 * @brief Main loop of the Nelder Mead algorithm
 *
 */
void opt::step()
{
    double f_1 = f(points[0], context);
    double f_n = f(points[points.size() - 2], context);
    double f_n1 = f(points[points.size() - 1], context);
    vector<double> centroid = get_centroid();

    // Step 1
    sort_by_opt_function();
#ifdef DEBUG
    print_points();
#endif
    // Step 2
    double f_r;
    vector<double> x_r(d, 0.);
    if (reflect(centroid, f_1, f_n, f_r, x_r))
    {
        return;
    }
    f_r = isnan(f_r) ? DBL_MAX : f_r;
#ifdef DEBUG
    cout << "Reflection point ";
    vector_println(x_r);
#endif

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
#ifdef DEBUG
        cout << "outside contraction point " << x_oc[0] << "\n";
#endif
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
#ifdef DEBUG
    cout << "f(x_r) = " << f_r << ", f(x_n+1) = " << f_n1 << "\n";
#endif

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

void opt::print_points()
{
    cout << "[\n";
    for (auto it = points.begin(); it != points.end(); it++)
    {
        cout << "  [";
        for (auto v = (*it).begin(); v != (*it).end(); v++)
        {
            cout << *v << ", ";
        }
        cout << "]\n";
    }
    cout << "]\n";
}

/**
 * @brief Constructs an optimization solution using the first point in `points`
 *
 * @return solution* A holder for the optimal vector at the current state of
 * `opt`, and the value of the objective function applied to this vector
 */
solution *opt::make_solution()
{
    return new solution(points.front(), f(points.front(), context));
}

/**
 * @brief Evaluates the objective function on all points and returns the results
 * in a vector.
 *
 * @return vector<double> [f(x_1), f(x_2), ..., f(x_n)]
 */
vector<double> opt::eval_all()
{
    vector<double> output(points.size(), 0.);
    for (int i = 0; i < output.size(); i++)
    {
        output[i] = f(points[i], context);
    }
    return output;
}

/**
 * @brief Termination condition for the Nelder-Mead algorithm, determined using
 * the change in sample standard deviation of the points in the polytope
 *
 * @return true iff the fractional change in sample standard deviation since the
 * last call to this function is less than `STDDEV_TOL`.
 */
bool opt::should_terminate()
{
    vector<double> vals = eval_all();
    int n = vals.size();
    double mean = 0;
    for (int i = 0; i < n; i++)
    {
        mean += vals[i] / n;
    }

    double stddev = 0;
    for (int i = 0; i < n; i++)
    {
        stddev += pow(vals[i] - mean, 2.) / n;
    }
    bool output = abs(stddev / last_stddev - 1) < STDDEV_TOL;
    output = output || stddev == last_stddev; // Handle 0/0
    last_stddev = stddev;
    return output;
}

solution *opt::solve()
{
    if (points.size() < d + 1)
    {
        throw invalid_argument("Insufficient initial points provided before calling `opt::solve`");
    }
    bool started = false;
    while (!should_terminate() || !started)
    {
        step();
        started = true;
    }
    return make_solution();
}

// Test functions
double fun(const vector<double> &x, void *ctx)
{
    return pow((x[0] - 3), 2) + pow((x[1] + 5), 2) + pow((x[2] - 12), 2);
}

typedef struct fixed_temp
{
    ModifiedRaoultModel m;
    double temp;
} * fixed_temp_t;

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

// int main()
// {
//     AntoineModel a1 = ANTOINE_METHANOL;
//     AntoineModel a2 = ANTOINE_WATER;
//     BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398, 363.15, T_unit::K);
//     ModifiedRaoultModel m(a1, a2, b, T_unit::K, 1, P_unit::bar);
//     double meth_bp = a1.tsat(1);
//     double water_bp = a2.tsat(1);

//     double query_temperature = 350;

//     // Solve for x given T
//     fixed_temp_t context = (fixed_temp_t)malloc(sizeof(fixed_temp));
//     context->m = m;
//     context->temp = query_temperature;

//     opt o(p_error, context, 1);
//     vector<vector<double>> initial_points = zip(
//         vector<vector<double>>{linspace(0, 1, 10)});
//     o.set_polytope(initial_points);

//     solution *s = o.solve();
//     cout << (s->x)[0] << '\n';
// }

int main()
{
    opt o(fun, NULL, 3);
    vector<vector<double>> initial_points = zip(vector<vector<double>>{
        linspace(-20, 20, 5),
        linspace(-20, 20, 5),
        linspace(-20, 20, 5),
    });

    o.set_polytope(initial_points);
    solution *s = o.solve();
    vector_println(s->x);
    return 0;
}