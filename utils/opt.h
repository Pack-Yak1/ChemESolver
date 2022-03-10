#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

#ifndef opt_h
#define opt_h

static const double ALPHA = 1, BETA = 2, GAMMA = 0.5, DELTA = 0.5;

typedef double (*opt_func_t)(const vector<double> &x);

class opt
{
private:
    opt_func_t f;
    unsigned int d;

    void sort_by_opt_function(vector<vector<double>> &points);
    vector<double> get_centroid(const vector<vector<double>> &points);
    double reflect(vector<vector<double>> &points, const vector<double> &centroid, double f_1, double f_n);

public:
    opt(opt_func_t f, unsigned int d);
};

#endif