#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>

#include "vector_ops.hpp"

using namespace std;

#ifndef NM_opt
#define NM_opt

const int NON_TERMINATING = 1e5;

const double STDDEV_TOL = 1e-30;

typedef double (*opt_func_t)(const vector<double> &x, void *context);

class solution
{
public:
    vector<double> x;
    double fx;

private:
    solution(vector<double> x, double fx);
    friend class opt;
};

class opt
{
private:
    opt_func_t f;
    size_t d;
    vector<vector<double>> points;
    vector<double> fx_cache;
    vector<double> centroid;
    vector<double> last_accepted;
    void *context;
    unsigned int num_points;
    double ALPHA;
    double BETA;
    double GAMMA;
    double DELTA;

    void accept(vector<double> point, double value);
    void sort_by_opt_function(bool accepted);
    void set_centroid(bool accepted);
    bool reflect(const vector<double> &centroid, double f_1, double f_n,
                 double &f_r, vector<double> &x_r);
    void expand(const vector<double> &centroid, const vector<double> &x_r,
                double f_r);
    vector<double> outside_contract(const vector<double> &centroid,
                                    const vector<double> &x_r);
    bool inside_contract(const vector<double> &centroid, double f_n1);
    void shrink(const vector<double> &x1);
    bool step(bool accepted);
    solution *make_solution();
    vector<double> eval_all();
    bool should_terminate();
    void set_polytope(vector<vector<double>> &points);
    solution *solve_helper();
    void print_points(bool display_fx = false);

public:
    /**
     * @brief Construct a new optimization problem.
     *
     * @param f The optimization function to be minimized
     * @param context Pointer to additional data needed for the optimization
     * function
     * @param d The dimensionality of the optimization problem
     */
    opt(opt_func_t f, void *context, unsigned int d);

    /**
     * @brief Setter method for context of optimization problem
     *
     * @param context The context which will be used by the optimization
     * function provided
     */
    void set_context(void *context);

    /**
     * @brief Run the Nelder-Mead algorithm on the provided optimization
     * problem. Overwrites the current polytope and starts the algorithm with
     * a random polytope consisting points where each coordinate is between min
     * and max. For better performance, min and max should span the possible
     * domain of the optimal solution
     * @param lb An estimate for the lower bounds of solutions.
     * @param ub An estimate for the upper bounds of solutions.
     *
     * @return solution* containing the found optimal point and the value of the
     * objective function at this point. Caller is responsible for deleting this
     * solution object.
     */
    solution *solve(const vector<double> &lb, const vector<double> &ub);
};

#endif