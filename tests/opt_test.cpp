#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "../utils/opt.h"
#include <chrono>

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
    const double MIN = -100;
    const double MAX = 100;
    const double NUM_ITERS = 10;
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(MIN, MAX);
    for (int num_dimensions = 1; num_dimensions < 5; num_dimensions++)
    {
        cout << "Test solving a " << num_dimensions << "-dimensional problem with random solution\n";
        double max_error = 0.;
        double total_time = 0.;
        double total_error = 0.;
        for (int test_iter = 0; test_iter < NUM_ITERS; test_iter++)
        {
            // Generate random answer set
            vector<double> answers(num_dimensions, 0.);
            for (int i = 0; i < num_dimensions; i++)
            {
                answers[i] = distr(eng);
            }
            cout << "Test iteration " << test_iter << " for " << num_dimensions << "-dimensional problem\n";
            cout << "True solution:     ";
            vector_println(answers);

            // Produce solution and record time
            opt o2(foo, &answers, num_dimensions);
            auto start = chrono::high_resolution_clock::now();
            solution *s2 = o2.auto_solve(MIN, MAX);
            auto stop = chrono::high_resolution_clock::now();
            total_time += chrono::duration<double>(stop - start).count();

            // Compare solution to answers
            cout << "Solution obtained: ";
            vector_println(s2->x);
            cout << "Value of objective function: " << s2->fx << "\n";
            vector<double> errs(num_dimensions, 0.);
            for (int i = 0; i < num_dimensions; i++)
            {
                errs[i] = abs(answers[i] - s2->x[i]);
            }
            cout << "Errors:            ";
            vector_println(errs);
            cout << '\n';

            // Update max error seen and total error
            double iter_max_error = *max_element(
                errs.begin(), errs.end());
            max_error = max(iter_max_error, max_error);
            total_error += vector_mean(errs);
        }
        cout << "Average time taken: " << total_time / NUM_ITERS << "\n";
        cout << "Highest error: " << max_error << "\n";
        cout << "Average error: " << total_error / NUM_ITERS << "\n\n";
    }
    return 0;
}