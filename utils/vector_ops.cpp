#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "vector_ops.h"

using namespace std;

vector<double> sum(vector<double> v1, vector<double> v2)
{
    if (v1.size() != v2.size())
    {
        throw DIMENSION_ERROR;
    }
    int d = v1.size();
    vector<double> output(d, 0.);
    for (auto dim = 0; dim < v1.size(); dim++)
    {
        output[dim] += v1[dim] + v2[dim];
    }
    return output;
}

vector<double> sum(const vector<vector<double>> &vecs)
{
    if (vecs.size() == 0)
    {
        return vector<double>();
    }
    int d = vecs.front().size();
    vector<double> output(d, 0.);
    for (auto it = vecs.begin(); it != vecs.end(); it++)
    {
        if (d != it->size())
        {
            throw DIMENSION_ERROR; // All vectors in vecs must have dimensions d
        }
        for (auto dim = 0; dim < it->size(); dim++)
        {
            output[dim] += (*it)[dim];
        }
    }
    return output;
}

/**
 * @brief Returns a vector of linearly spaced points
 *
 * @param start The lower bound for the linear space, inclusive
 * @param end The upper bound for the linear space, inclusive
 * @param points The number of points in the resulting linear space
 * @return vector<double> A linear space of doubles starting at `start`, ending
 * at `end`, with `points` points.
 */
vector<double> linspace(double start, double end, double points)
{
    if (points <= 0)
    {
        return vector<double>();
    }
    vector<double> output(points, 0.);
    double step_size = (end - start) / (points - 1);
    output[0] = start;
    for (int i = 1; i < output.size(); i++)
    {
        output[i] = output[i - 1] + step_size;
    }
    return output;
}

/**
 * @brief Similar to numpy meshgrid. Given n linspaces, return a vector of n-
 * dimensional vectors v_i spanning the cartesian product of linspaces
 *
 * @param linspaces A vector of linspaces which the resulting grid must span
 * @return vector<vector<double>> The resulting meshgrid
 */
vector<vector<double>> zip(vector<vector<double>> linspaces)
{
    int num_dimensions = linspaces.size();
    if (num_dimensions == 0)
    {
        return vector<vector<double>>();
    }

    vector<vector<double>> prev;
    vector<vector<double>> next;
    for (int i = 0; i < linspaces[0].size(); i++)
    {
        next.push_back(vector<double>{linspaces[0][i]});
    }

    for (int i = 1; i < num_dimensions; i++)
    {
        prev = next;
        next.clear();
        // From each vector v in prev, add all vectors starting with v and
        // ending with a value from linspaces[i].
        for (int j = 0; j < prev.size(); j++)
        {
            for (int k = 0; k < linspaces[i].size(); k++)
            {
                vector<double> to_add = prev[j];
                to_add.emplace_back(linspaces[i][k]);
                next.push_back(to_add);
            }
        }
    }
    return next;
}

// int main()
// {
//     vector<double> tmp = linspace(1., 10., 10.);
//     vector<double> tmp2 = linspace(11., 20., 10.);
//     vector<vector<double>> gg{tmp, tmp2};
//     vector<vector<double>> test = zip(gg);

//     for (int i = 0; i < test.size(); i++)
//     {
//         for (int j = 0; j < test[0].size(); j++)
//         {
//             cout << test[i][j] << ", ";
//         }
//         cout << "\n";
//     }
// }
