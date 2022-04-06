#include "opt.hpp"
#include "vector_ops.hpp"

#include <float.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

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
    this->points = vector<vector<double>>();
    this->num_points = 0;
    this->fx_cache = vector<double>(d, 0.);
    this->centroid = vector<double>(d, 0.);
    this->last_accepted = vector<double>(d, 0.);

    this->ALPHA = 1;
#ifdef ANMS
    this->BETA = d >= 2 ? 1 + 2 / d : 2.;
    this->GAMMA = d >= 2 ? 0.75 - 0.5 / d : 0.5;
    this->DELTA = d >= 2 ? 1. - 1. / d : 0.5;
#else
    this->BETA = 2.;
    this->GAMMA = 0.5;
    this->DELTA = 0.5;
#endif
}

void opt::accept(vector<double> point, double value)
{
    points.back() = point;
    if (isinf(value) || isnan(value))
    {
        value = DBL_MAX;
    }
    fx_cache.back() = value;
    last_accepted = point;
}

/**
 * @brief Set the points defining the simplex used for the Nelder-Mead
 * algorithm
 *
 * @param points A vector of vectors [v1, v2, ..., vn] such that all vi are
 * vectors of size `d`.
 */
void opt::set_polytope(vector<vector<double>> &points)
{
    this->points = points;
    this->num_points = points.size();
}

void opt::set_context(void *context) { this->context = context; }

// void swap(vector<double> &arr, int i, int j)
// {
//     double tmp = arr[i];
//     arr[i] = arr[j];
//     arr[j] = tmp;
// }

template <typename T>
void swap(vector<T> &arr, int i, int j)
{
    T tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

/**
 * @brief Reorders vectors v1 and v2 by in the order specified by indices.
 * Assumes indices contains exactly one of each integer from [0..n-1] where n
 * is the size of all arguments passed to this function.
 *
 * Destructively modifies `indices`. Values of `indices` are unspecified after
 * this function is called.
 *
 */
template <typename T1, typename T2>
void apply_indices(vector<double> &indices, vector<T1> &v1, vector<T2> &v2)
{
    if (!(indices.size() == v1.size() && v1.size() == v2.size()))
    {
        throw invalid_argument(
            "Attempted to apply indices to vectors of different size");
    }
    size_t n = indices.size();
    size_t tmp;
    for (size_t origin = 0; origin < n; origin++)
    {
        size_t i = origin;
        while (origin != indices[i] && indices[i] != -1)
        {
            swap(v1, indices[i], i);
            swap(v2, indices[i], i);
            tmp = indices[i];
            indices[i] = -1;
            i = tmp;
        }
        indices[i] = -1;
    }
}

/**
 * @brief Sorts `points` in place. After sorting, `points` = `[x1, ..., xn]`
 * such that f(x1) <= ... <= f(xn).
 *
 * @param accepted True if this function is called in an iteration of the
 * Nelder-Mead algorithm after a non-shrink step which invoked `accept`.
 */
#define noop (void)0
void opt::sort_by_opt_function(bool accepted)
{
    if (accepted)
    {
        // Manual search without std::upper_bound to handle inf/nan.
        // Insertion sort on last element for O(n) sort
        double to_insert = fx_cache.back();
        if (isnan(to_insert) || isinf(to_insert))
        {
            return;
        }
        size_t idx = 0;
        while (idx < num_points - 1 && fx_cache[idx] <= to_insert && !isnan(fx_cache[idx]) && !isnan(fx_cache[idx]))
        {
            idx++;
        }
        auto cache_start = fx_cache.begin() + idx;
        std::rotate(cache_start, fx_cache.end() - 1, fx_cache.end());

        auto points_start = points.begin() + idx;
        std::rotate(points_start, points.end() - 1, points.end());
    }
    else
    {
        // Comparator to sort by cached value in fx_cache
        auto comparator = [this](size_t left, size_t right) -> bool
        {
            double cmp;
            double f1 = fx_cache[left];
            double f2 = fx_cache[right];
            if (isnan(f1) || isinf(f1))
            {
                cmp = 1;
            }
            else if (isnan(f2) || isinf(f2))
            {
                cmp = -1;
            }
            else
            {
                cmp = f1 - f2;
            }
            if (cmp == 0)
            {
                // Consistent tie break needed for Nelder-Meads
                for (size_t i = 0; i < d; i++)
                {
                    double entry_cmp = points[left][i] - points[right][i];
                    if (entry_cmp != 0)
                    {
                        return entry_cmp < 0;
                    }
                }
            }
            return cmp < 0;
        };
        // Argsort and store indices
        vector<double> indices(num_points);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), comparator);

        apply_indices(indices, points, fx_cache);
    }
}

/**
 * @brief Sets the centroid of all but the last point in `points`. If `accepted`
 * is true, assumes `points` is sorted by `sort by opt function` and that
 * `last_accepted` is set.
 *
 * @param accepted True if this function is called in an iteration of the
 * Nelder-Mead algorithm after a non-shrink step which invoked `accept`.
 *
 */
void opt::set_centroid(bool accepted)
{
    double denom = (double)num_points - 1.;
    if (accepted)
    {
        sum_in_place(centroid, points.back(), 1., -1. / denom);
        sum_in_place(centroid, last_accepted, 1., 1. / denom);
        return;
    }
    centroid = vector<double>(d, 0.);
    for (unsigned int i = 0; i < num_points - 1; i++)
    {
        sum_in_place(centroid, points[i], 1., 1. / denom);
#ifdef DEBUG
        cout << "iter " << i << " ended. output = ";
        vector_println(centroid);
#endif
    }
#ifdef DEBUG
    cout << "Centroid calculated: ";
    vector_println(centroid);
#endif
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
bool opt::reflect(const vector<double> &centroid, double f_1, double f_n,
                  double &f_r, vector<double> &x_r)
{
    vector<double> last = points.back();
    sum_in_place(x_r, multiply(1 + ALPHA, centroid));
    sum_in_place(x_r, multiply(-ALPHA, last));
#ifdef DEBUG
    cout << "Reflection point calculated, points were";
    print_points(true);
    cout << "Reflection point is";
    vector_println(x_r);
#endif
    f_r = f(x_r, context);
    if (f_1 <= f_r && f_r < f_n)
    {
#ifdef DEBUG
        cout << "Accepted reflection point\n";
#endif
        accept(x_r, f_r);
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
void opt::expand(const vector<double> &centroid, const vector<double> &x_r,
                 double f_r)
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
#ifdef DEBUG
        cout << "accepted expansion point\n";
#endif
        accept(x_e, f_e);
    }
    else
    {
#ifdef DEBUG
        cout << "accepted reflection point\n";
#endif
        accept(x_r, f_r);
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
vector<double> opt::outside_contract(const vector<double> &centroid,
                                     const vector<double> &x_r)
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
    sum_in_place(x_ic, multiply(1 + GAMMA, centroid));
    sum_in_place(x_ic, multiply(-GAMMA, points.back()));
#ifdef DEBUG
    cout << "inside contraction point ";
    vector_println(x_ic);
#endif
    double f_ic = f(x_ic, context);
    if (f_ic < f_n1)
    {
#ifdef DEBUG
        cout << "accepted inside contraction point\n";
#endif
        accept(x_ic, f_ic);
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
    for (unsigned int i = 1; i < num_points; i++)
    {
        multiply_in_place(DELTA, points[i]);
        sum_in_place(points[i], multiply(1 - DELTA, x_1));
        double result = f(points[i], context);
        fx_cache[i] = isinf(result) || isnan(result) ? DBL_MAX : result;
        // fx_cache[i] = f(points[i], context);
#ifdef DEBUG
        cout << "Point " << i << " was updated to ";
        vector_println(points[i]);
#endif
    }
}

/**
 * @brief Main loop of the Nelder Mead algorithm
 *
 * @return true iff this iteration accepted a single point (i.e. did not use a
 * shrink operation)
 *
 */
bool opt::step(bool accepted)
{
    // print_points(true);
#ifdef DEBUG
    cout << "Step started with points: \n";
    print_points(true);
#endif
    // Step 1
    sort_by_opt_function(accepted);
#ifdef DEBUG
    cout << "Sorted points: \n";
    print_points(true);
#endif
    double f_1 = fx_cache.front();
    double f_n = fx_cache[num_points - 2];
    double f_n1 = fx_cache.back();
    set_centroid(accepted);

    // Step 2
    double f_r;
    vector<double> x_r(d, 0.);
    if (reflect(centroid, f_1, f_n, f_r, x_r))
    {
        return true;
    }

    // Step 3
    if (f_r < f_1)
    {
        expand(centroid, x_r, f_r);
        return true;
    }

    // Step 4
    vector<double> x_oc;
    if (f_n <= f_r && f_r < f_n1)
    {
        x_oc = outside_contract(centroid, x_r);
#ifdef DEBUG
        cout << "outside contraction point ";
        vector_println(x_oc);
#endif
        double f_oc = f(x_oc, context);
        if (f_oc <= f_r)
        {
#ifdef DEBUG
            cout << "accepted outside contraction point \n";
#endif
            accept(x_oc, f_oc);
            return true;
        }
        else
        {
            // Skip to step 6
            shrink(points[0]);
            return false;
        }
    }

    // Step 5
    if (f_r >= f_n1)
    {
        if (inside_contract(centroid, f_n1))
        {
            return true;
        }
    }

    // Step 6
    shrink(points[0]);
    return false;
}

/**
 * @brief Pretty print the vectors in `points`
 *
 * @param display_fx If true, prints the value of the objective function next
 * to each point printed.
 */
void opt::print_points(bool display_fx)
{
    cout << "[\n";
    for (auto it = points.begin(); it != points.end(); it++)
    {
        cout << "  [";
        for (auto v = (*it).begin(); v != (*it).end(); v++)
        {
            cout << *v << ", ";
        }
        cout << "]";
        if (display_fx)
        {
            cout << "     f(x) = " << fx_cache[(size_t)(it - points.begin())] << ".";
        }
        cout << "\n";
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
    return new solution(points.front(), fx_cache.front());
}

/**
 * @brief Evaluates the objective function on all points and returns the results
 * in a vector.
 *
 * @return vector<double> [f(x_1), f(x_2), ..., f(x_n)]
 */
vector<double> opt::eval_all()
{
    vector<double> output(num_points, 0.);
    for (size_t i = 0; i < output.size(); i++)
    {
        if (points[i].size() != d)
        {
            throw invalid_argument(
                "Dimension of vector do not match specified dimension in "
                "`eval_all`. Ensure that all points provided to `set_polytope`"
                " match the dimension of the problem.");
        }
        double result = f(points[i], context);
        output[i] = isinf(result) || isnan(result) ? DBL_MAX : result;
        // output[i] = f(points[i], context);
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
    // int n = vals.size();
    double mean = 0;
    for (unsigned int i = 0; i < num_points; i++)
    {
        mean += fx_cache[i] / num_points;
    }

    double stddev = 0;
    for (unsigned int i = 0; i < num_points; i++)
    {
        stddev += pow(fx_cache[i] - mean, 2.) / num_points;
    }
    // cout << "stddev = " << stddev << '\n';
    bool output = stddev < STDDEV_TOL;
    return output;
}

/**
 * @brief Run the Nelder-Mead algorithm on the provided optimization
 * problem. Requires that user has set the initial polytope
 *
 * @return solution* containing the found optimal point and the value of the
 * objective function at this point.
 */
solution *opt::solve_helper()
{
    if (num_points < d + 1)
    {
        throw invalid_argument(
            "Insufficient initial points provided before calling `opt::solve`");
    }
    int num_iters = 0;
    bool accepted = false;
    while (!should_terminate())
    {
        accepted = step(accepted);
        num_iters++;
#ifdef DEBUG
        cout << "\n=============="
             << "ITERATION NUMBER " << num_iters << "=================\n";
#endif
        if (num_iters > NON_TERMINATING)
        {
            cerr
                << "Algorithm failed to converge in " << NON_TERMINATING
                << " steps. Returning best point found at current iteration.\n";
            return make_solution();
        }
    }
    return make_solution();
}

solution *opt::solve(const vector<double> &lb, const vector<double> &ub)
{
    vector<vector<double>> initial_points;
    double num_vertices = d + 1;
    initial_points.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++)
    {
        vector<double> to_add(lb);
        if (i < d)
        {
            to_add[i] = ub[i];
        }
        initial_points.emplace_back(to_add);
    }
    set_polytope(initial_points);
#ifdef DEBUG
    cout << "solver started with points:\n";
    print_points();
#endif
    this->fx_cache = eval_all();
    set_centroid(false);
    return solve_helper();
}
