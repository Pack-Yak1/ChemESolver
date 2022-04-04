#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

#ifndef vector_ops
#define vector_ops

/**
 * @brief Returns the mean of v
 */
double vector_mean(const vector<double> &v);

/**
 * @brief Multiplies `v` in place by `c`.
 */
void multiply_in_place(double c, vector<double> &v);

/**
 * @brief Returns `v` multiplied by scalar constant `c`.
 */
vector<double> multiply(double c, const vector<double> &v);

/**
 * @brief Stores the vector sum of `c1 * v1` and `c2 * v2` in `v1`.
 */
void sum_in_place(vector<double> &v1, const vector<double> &v2, double c1 = 1.,
                  double c2 = 1.);

vector<double> sum(const vector<double> &v1, const vector<double> &v2);

vector<double> sum(const vector<vector<double>> &vecs);

/**
 * @brief Returns a vector of linearly spaced points
 *
 * @param start The lower bound for the linear space, inclusive
 * @param end The upper bound for the linear space, inclusive
 * @param points The number of points in the resulting linear space
 * @return vector<double> A linear space of doubles starting at `start`, ending
 * at `end`, with `points` points.
 */
vector<double> linspace(double start, double end, double points);

/**
 * @brief Similar to numpy meshgrid. Given n linspaces, return a vector of n-
 * dimensional vectors v_i spanning the cartesian product of linspaces
 *
 * @param linspaces A vector of linspaces which the resulting grid must span
 * @return vector<vector<double>> The resulting meshgrid
 */
vector<vector<double>> zip(const vector<vector<double>> &linspaces);

void vector_print(const vector<double> &v);

void vector_println(const vector<double> &v);

void matrix_print(const vector<vector<double>> &m);

#endif