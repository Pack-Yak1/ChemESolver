#include "opt.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "vector_ops.h"

using namespace std;

solution::solution(vector<double> x, double fx) {
    this->fx = fx;
    this->x = x;
}

opt::opt(opt_func_t f, void *context, unsigned int d) {
    this->f = f;
    this->d = d;
    this->context = context;
    this->last_stddev = 0;
    this->points = vector<vector<double>>();
    this->num_points = 0;
}

/**
 * @brief Set the points defining the simplex used for the Nelder-Mead
 * algorithm
 *
 * @param points A vector of vectors [v1, v2, ..., vn] such that all vi are
 * vectors of size `d`.
 */
void opt::set_polytope(vector<vector<double>> &points) {
    this->points = points;
    this->num_points = points.size();
}

void opt::set_context(void *context) { this->context = context; }

/**
 * @brief Sorts `points` in place. After sorting, `points` = `[x1, ..., xn]`
 * such that f(x1) <= ... <= f(xn).
 */
void opt::sort_by_opt_function() {
    auto comparator = [&](vector<double> p1, vector<double> p2) -> bool {
        // Check if dimensions correct
        if (p1.size() != d || p2.size() != d) {
            // `sort_by_opt_function` given `points` argument containing vector
            // of incorrect dimensions.
            throw invalid_argument(
                "Dimensions of vectors do not match in `sort_by_opt_function`. "
                "Ensure that all points provided to `set_polytope` match the "
                "dimension of the problem.");
        }
        double f1 = f(p1, context);
        double f2 = f(p2, context);
        double cmp;
        if (isnan(f1)) {
            cmp = 1;
        } else if (isnan(f2)) {
            cmp = -1;
        } else {
            cmp = f(p1, context) - f(p2, context);
        }
        if (cmp == 0) {
            // Consistent tie break needed for Nelder-Meads
            for (int i = 0; i < d; i++) {
                double entry_cmp = p1[i] - p2[i];
                if (entry_cmp != 0) {
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
vector<double> opt::get_centroid() {
    vector<double> output(d, 0.);
    for (int i = 0; i < num_points - 1; i++) {
        sum_in_place(output, points[i]);
#ifdef DEBUG
        cout << "iter " << i << " ended. output = ";
        vector_println(output);
#endif
    }
    multiply_in_place(1 / ((double)num_points - 1), output);
#ifdef DEBUG
    cout << "Centroid calculated: ";
    vector_println(output);
#endif
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
bool opt::reflect(const vector<double> &centroid, double f_1, double f_n,
                  double &f_r, vector<double> &x_r) {
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
    if (f_1 <= f_r && f_r < f_n) {
#ifdef DEBUG
        cout << "Accepted reflection point\n";
#endif
        points.back() = x_r;
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
                 double f_r) {
    vector<double> x_e(d, 0.);
    sum_in_place(x_e, multiply(1 - BETA, centroid));
    sum_in_place(x_e, multiply(BETA, x_r));
    double f_e = f(x_e, context);
#ifdef DEBUG
    cout << "Expansion point ";
    vector_println(x_e);
#endif
    if (f_e < f_r) {
#ifdef DEBUG
        cout << "accepted expansion point\n";
#endif
        points.back() = x_e;
    } else {
#ifdef DEBUG
        cout << "accepted reflection point\n";
#endif
        points.back() = x_r;
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
                                     const vector<double> &x_r) {
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
bool opt::inside_contract(const vector<double> &centroid, double f_n1) {
    vector<double> x_ic(d, 0.);
    sum_in_place(x_ic, multiply(1 + GAMMA, centroid));
    sum_in_place(x_ic, multiply(-GAMMA, points.back()));
#ifdef DEBUG
    cout << "inside contraction point ";
    vector_println(x_ic);
#endif
    double f_ic = f(x_ic, context);
    if (f_ic < f_n1) {
#ifdef DEBUG
        cout << "accepted inside contraction point\n";
#endif
        points.back() = x_ic;
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
void opt::shrink(const vector<double> &x_1) {
    for (int i = 1; i < num_points; i++) {
        multiply_in_place(DELTA, points[i]);
        sum_in_place(points[i], multiply(1 - DELTA, x_1));
#ifdef DEBUG
        cout << "Point " << i << " was updated to ";
        vector_println(points[i]);
#endif
    }
}

/**
 * @brief Main loop of the Nelder Mead algorithm
 *
 */
void opt::step() {
#ifdef DEBUG
    cout << "Step started with points: \n";
    print_points(true);
#endif
    // Step 1
    sort_by_opt_function();
#ifdef DEBUG
    cout << "Sorted points: \n";
    print_points(true);
#endif
    double f_1 = f(points.front(), context);
    double f_n = f(points[num_points - 2], context);
    double f_n1 = f(points.back(), context);
    vector<double> centroid = get_centroid();

    // Step 2
    double f_r;
    vector<double> x_r(d, 0.);
    if (reflect(centroid, f_1, f_n, f_r, x_r)) {
        return;
    }
    f_r = isnan(f_r) ? DBL_MAX : f_r;

    // Step 3
    if (f_r < f_1) {
        expand(centroid, x_r, f_r);
        return;
    }

    // Step 4
    vector<double> x_oc;
    if (f_n <= f_r && f_r < f_n1) {
        x_oc = outside_contract(centroid, x_r);
#ifdef DEBUG
        cout << "outside contraction point ";
        vector_println(x_oc);
#endif
        double f_oc = f(x_oc, context);
        f_oc = isnan(f_oc) ? DBL_MAX : f_oc;
        if (f_oc <= f_r) {
#ifdef DEBUG
            cout << "accepted outside contraction point \n";
#endif
            points.back() = x_oc;
            return;
        } else {
            // Skip to step 6
            shrink(points[0]);
            return;
        }
    }

    // Step 5
    if (f_r >= f_n1) {
        if (inside_contract(centroid, f_n1)) {
            return;
        }
    }

    // Step 6
    shrink(points[0]);
}

/**
 * @brief Pretty print the vectors in `points`
 *
 * @param display_fx If true, prints the value of the objective function next
 * to each point printed.
 */
void opt::print_points(bool display_fx) {
    cout << "[\n";
    for (auto it = points.begin(); it != points.end(); it++) {
        cout << "  [";
        for (auto v = (*it).begin(); v != (*it).end(); v++) {
            cout << *v << ", ";
        }
        cout << "]";
        if (display_fx) {
            cout << "     " << f(*it, context);
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
solution *opt::make_solution() {
    sort_by_opt_function();
    return new solution(points.front(), f(points.front(), context));
}

/**
 * @brief Evaluates the objective function on all points and returns the results
 * in a vector.
 *
 * @return vector<double> [f(x_1), f(x_2), ..., f(x_n)]
 */
vector<double> opt::eval_all() {
    vector<double> output(num_points, 0.);
    for (int i = 0; i < output.size(); i++) {
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
bool opt::should_terminate() {
    vector<double> vals = eval_all();
    int n = vals.size();
    double mean = 0;
    for (int i = 0; i < n; i++) {
        mean += vals[i] / n;
    }

    double stddev = 0;
    for (int i = 0; i < n; i++) {
        stddev += pow(vals[i] - mean, 2.) / n;
    }
    bool output = abs(stddev / last_stddev - 1) < STDDEV_TOL;
    output = output || stddev == last_stddev;  // Handle 0/0
    last_stddev = stddev;
    return output;
}

/**
 * @brief Run the Nelder-Mead algorithm on the provided optimization
 * problem. Requires that user has set the initial polytope
 *
 * @return solution* containing the found optimal point and the value of the
 * objective function at this point.
 */
solution *opt::solve_helper() {
    if (num_points < d + 1) {
        throw invalid_argument(
            "Insufficient initial points provided before calling `opt::solve`");
    }
    int num_iters = 0;
    bool started = false;
    while (!should_terminate() || !started || num_iters < MIN_ITERS) {
        step();
        started = true;
        num_iters++;
        if (num_iters > NON_TERMINATING) {
            // throw runtime_error("Algorithm failed to converge\n");
            cout
                << "Algorithm failed to converge in " << NON_TERMINATING
                << " steps. Returning best point found at current iteration.\n";
            return make_solution();
        }
    }
    return make_solution();
}

solution *opt::solve(const vector<double> &lb, const vector<double> &ub) {
    vector<vector<double>> initial_points;
    double num_vertices = d + 1;
    initial_points.reserve(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        vector<double> to_add(lb);
        if (i < d) {
            to_add[i] = ub[i];
        }
        initial_points.emplace_back(to_add);
    }
    set_polytope(initial_points);
    return solve_helper();
}
