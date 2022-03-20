#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "vector_ops.h"

using namespace std;

#ifndef opt_h
#define opt_h

#define STDDEV_TOL 1e-10

static const double ALPHA = 0.5, BETA = 2, GAMMA = .5, DELTA = 0.5;

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
    unsigned int d;
    vector<vector<double>> points;
    void *context;
    double last_stddev;

    void sort_by_opt_function();
    vector<double> get_centroid();
    bool reflect(const vector<double> &centroid, double f_1, double f_n, double &f_r, vector<double> &x_r);
    void expand(const vector<double> &centroid, const vector<double> &x_r, double f_r);
    vector<double> outside_contract(const vector<double> &centroid, const vector<double> &x_r);
    bool inside_contract(const vector<double> &centroid, double f_n1);
    void shrink(const vector<double> &x1);
    void step();
    solution *make_solution();
    vector<double> eval_all();
    bool should_terminate();

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
     * @brief Set the points defining the simplex used for the Nelder-Mead
     * algorithm
     *
     * @param points A vector of vectors [v1, v2, ..., vn] such that all vi are
     * vectors of size `d`.
     */
    void set_polytope(vector<vector<double>> &points);

    /**
     * @brief Pretty print the vectors in `points`
     */
    void print_points();

    /**
     * @brief Setter method for context of optimization problem
     *
     * @param context The context which will be used by the optimization
     * function provided
     */
    void set_context(void *context);

    /**
     * @brief Run the Nelder-Mead algorithm on the provided optimization
     * problem.
     *
     * @return solution* containing the optimal point and the value of the
     * objective function at this point.
     */
    solution *solve();
};

#endif