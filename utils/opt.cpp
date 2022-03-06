#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

typedef double (*opt_func_t)(vector<double> &x);

void sort_by_opt_function(opt_func_t f, vector<vector<double>> &points, int d)
{
    auto comparator = [&](vector<double> p1, vector<double> p2) -> bool
    {
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

double fun(vector<double> &x)
{
    return -pow((x[0] - 3), 2);
}

int main()
{
    vector<vector<double>> v;
    for (int i = 0; i < 10; i++)
    {
        v.emplace_back(vector<double>(1, i));
    }
    sort_by_opt_function(&fun, v, 1);
    for (int i = 0; i < 10; i++)
    {
        cout << v[i][0] << '\n';
    }
    return 0;
}