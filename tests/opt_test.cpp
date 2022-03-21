#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "../utils/opt.h"

double fun(const vector<double> &x, void *ctx)
{
    return pow((x[0] - 3), 2) + pow((x[1] + 5), 2) + pow((x[2] - 12), 2);
}

// Requires: f_data points to a vector of doubles with the same size as x.
double foo(const vector<double> &x, void *f_data)
{
    vector<double> answers = *(vector<double> *)f_data;
    double output = 0;
    int n = x.size();
    for (int i = 0; i < n; i++)
    {
        output += pow(x[i] - answers[i], 2.);
    }
    return output;
}

int main()
{
    // Test 1
    cout << "Test solving a 3-dimensional problem with solution [3, -5, 12]\n";
    opt o(fun, NULL, 3);
    vector<vector<double>> initial_points = zip(vector<vector<double>>{
        linspace(-20, 20, 3),
        linspace(-20, 20, 3),
        linspace(-20, 20, 3),
    });

    o.set_polytope(initial_points);
    solution *s = o.solve();
    cout << "Solution obtained: ";
    vector_println(s->x);
    cout << "Value of objective function: " << s->fx << "\n\n";

    // Test 2
    for (int NUM_DIMENSIONS = 1; NUM_DIMENSIONS < 10; NUM_DIMENSIONS++)
    {
        const double MIN = -100;
        const double MAX = 100;
        cout << "Test solving a " << NUM_DIMENSIONS << "-dimensional problem with random solution\n";

        random_device rd;
        default_random_engine eng(rd());
        uniform_real_distribution<double> distr(MIN, MAX);

        vector<double> answers(NUM_DIMENSIONS, 0.);
        for (int i = 0; i < NUM_DIMENSIONS; i++)
        {
            answers[i] = distr(eng);
        }

        cout << "True solution: ";
        vector_println(answers);
        opt o2(foo, &answers, NUM_DIMENSIONS);
        solution *s2 = o2.auto_solve(MIN, MAX);
        cout << "Solution obtained: ";
        vector_println(s2->x);
        cout << "Value of objective function: " << s2->fx << "\n";
        vector<double> percent_errs(NUM_DIMENSIONS, 0.);
        for (int i = 0; i < NUM_DIMENSIONS; i++)
        {
            percent_errs[i] = abs(answers[i] / s2->x[i] - 1) * 100;
        }
        cout << "Percent errors: ";
        vector_println(percent_errs);
        cout << '\n';
        cout << "Highest error: " << *max_element(percent_errs.begin(), percent_errs.end()) << "%\n\n";
    }

    return 0;
}