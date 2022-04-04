#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "vector_ops.hpp"

using namespace std;

double vector_mean(const vector<double> &v)
{
    double output = 0.;
    for (auto i = v.begin(); i != v.end(); i++)
    {
        output += *i;
    }
    return output / v.size();
}

void multiply_in_place(double c, vector<double> &v)
{
    int d = v.size();
    for (vector<double>::iterator i = v.begin(); i != v.end(); i++)
    {
        *i *= c;
    }
}

vector<double> multiply(double c, const vector<double> &v)
{
    int d = v.size();
    vector<double> output = v;
    for (vector<double>::iterator i = output.begin(); i != output.end(); i++)
    {
        *i *= c;
    }
    return output;
}

void sum_in_place(vector<double> &v1, const vector<double> &v2)
{
    if (v1.size() != v2.size())
    {
        throw invalid_argument("Dimensions of vectors do not match in `sum_in_place`.");
    }
    int d = v1.size();
    for (int dim = 0; dim < v1.size(); dim++)
    {
        v1[dim] += v2[dim];
    }
}

vector<double> sum(const vector<double> &v1, const vector<double> &v2)
{
    if (v1.size() != v2.size())
    {
        throw invalid_argument("Dimensions of vectors do not match in `sum`.");
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
            throw invalid_argument("Dimensions of vectors do not match in `sum`.");
        }
        for (auto dim = 0; dim < it->size(); dim++)
        {
            output[dim] += (*it)[dim];
        }
    }
    return output;
}

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

vector<vector<double>> zip(const vector<vector<double>> &linspaces)
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

void vector_print(const vector<double> &v)
{
    cout << '[';
    int len = v.size();
    for (int i = 0; i < len; i++)
    {
        cout << v[i];
        if (i < len - 1)
        {
            cout << ", ";
        }
    }
    cout << ']';
}

void matrix_print(const vector<vector<double>> &m)
{
    cout << '[';
    int len = m.size();
    for (int i = 0; i < len; i++)
    {
        cout << "  ";
        vector_print(m[i]);
        if (i < len - 1)
        {
            cout << ",\n";
        }
    }
    cout << "\n]\n";
}

void vector_println(const vector<double> &v)
{
    vector_print(v);
    cout << '\n';
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
