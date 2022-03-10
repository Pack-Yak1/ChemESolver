#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "opt.h"

using namespace std;

opt::opt(opt_func_t f, unsigned int d)
{
    this->f = f;
    this->d = d;
}

/**
 * @brief Sorts `points` in place. After sorting, `points` = `[x1, ..., xn]`
 * such that f(x1) <= ... <= f(xn).
 *
 * @param points The points to sort over
 */
void opt::sort_by_opt_function(vector<vector<double>> &points)
{
    auto comparator = [&](vector<double> p1, vector<double> p2) -> bool
    {
        // Check if dimensions correct
        if (p1.size() != d || p2.size() != d)
        {
            throw "`sort_by_opt_function` given `points` argument containing vector of incorrect dimensions.";
        }
        double cmp = f(p1) - f(p2);
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

vector<double> opt::get_centroid(const vector<vector<double>> &points)
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
 * @return The value of `f(x_r)` where `x_r` is the reflection point
 */
double opt::reflect(vector<vector<double>> &points, const vector<double> &centroid, double f_1, double f_n)
{
    vector<double> x_r(d, 0.);
    vector<double> last = points[points.size() - 1];
    for (int dim = 0; dim < d; dim++)
    {
        x_r[dim] = centroid[dim] * (1 + ALPHA) - last[dim];
    }
    double f_r = f(x_r);
    if (f_r >= f_1 && f_r < f_n)
    {
        points[points.size() - 1] = x_r;
    }
    return f_r;
}

// Test function
double fun(const vector<double> &x)
{
    return -pow((x[0] - 3), 2);
}

int main()
{
    opt o(fun, 1);
    vector<vector<double>> v;
    // for (int i = 0; i < 10; i++)
    // {
    //     v.emplace_back(vector<double>(1, i));
    // }
    // // o.sort_by_opt_function(v);
    // for (int i = 0; i < 10; i++)
    // {
    //     cout << v[i][0] << '\n';
    // }
    return 0;
}