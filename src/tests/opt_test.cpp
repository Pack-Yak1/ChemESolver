#include "opt.hpp"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <string.h>
#include <vector>

const int MIN_DIMENSIONS = 1;
const int MAX_DIMENSIONS = 101;
const double MIN = -100;
const double MAX = 100;
const double NUM_ITERS = 10;
const bool WRITE_TO_OUTPUT = true;

string STATS_FILENAME = "opt_stats.csv";

// Requires: f_data points to a vector of doubles with the same size as x.
double
foo(const vector<double> &x, void *f_data)
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

void n_dimensional_test_suite(int num_dimensions, int NUM_ITERS, double MIN,
                              double MAX, default_random_engine eng,
                              uniform_real_distribution<double> distr,
                              ostream &output_file)
{
    cout << "Test solving a " << num_dimensions
         << "-dimensional problem with random solution\n";
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
#ifdef VERBOSE
        cout << "Test iteration " << test_iter << " for " << num_dimensions
             << "-dimensional problem\n";
        cout << "True solution:     ";
        vector_println(answers);
#endif
        // Produce solution and record time
        opt o2(foo, &answers, num_dimensions);
        auto start = chrono::high_resolution_clock::now();
        solution *s2 = o2.solve(vector<double>(num_dimensions, MIN),
                                vector<double>(num_dimensions, MAX));
        auto stop = chrono::high_resolution_clock::now();
        total_time += chrono::duration<double>(stop - start).count();

        // Compare solution to answers
        vector<double> errs(num_dimensions, 0.);
        for (int i = 0; i < num_dimensions; i++)
        {
            errs[i] = abs(answers[i] / s2->x[i] - 1) * 100;
        }
#ifdef VERBOSE
        cout << "Solution obtained: ";
        vector_println(s2->x);
        cout << "Value of objective function: " << s2->fx << "\n";
        cout << "Errors:            ";
        vector_println(errs);
        cout << '\n';
#endif
        // Update max error seen and total error
        double iter_max_error = *max_element(errs.begin(), errs.end());
        max_error = max(iter_max_error, max_error);
        total_error += vector_mean(errs);
        delete s2;
    }
    cout << "Average time taken: " << total_time / NUM_ITERS << "\n";
    cout << "Highest error: " << max_error << "%\n";
    cout << "Average error: " << total_error / NUM_ITERS << "%\n\n";
    if (WRITE_TO_OUTPUT)
    {
        output_file << num_dimensions << "," << total_time / NUM_ITERS << "," << max_error << "," << total_error / NUM_ITERS << "\n";
    }
}

int main()
{
    // Test 1
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(MIN, MAX);
    ofstream output_file;
    if (WRITE_TO_OUTPUT)
    {
        output_file.open(STATS_FILENAME);
    }
    for (int num_dimensions = MIN_DIMENSIONS; num_dimensions < MAX_DIMENSIONS; num_dimensions++)
    {
        n_dimensional_test_suite(num_dimensions, NUM_ITERS, MIN, MAX, eng,
                                 distr, output_file);
    }
    if (WRITE_TO_OUTPUT)
    {
        output_file.close();
    }
    return 0;
}