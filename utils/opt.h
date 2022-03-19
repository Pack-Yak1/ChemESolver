#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

#ifndef opt_h
#define opt_h

#define DIMENSION_ERROR 1
#define DEGENERATE_POINTS 2

static const double ALPHA = 0.1, BETA = 1.1, GAMMA = .1, DELTA = 0.1;

typedef double (*opt_func_t)(const vector<double> &x, void *context);

class solution
{
public:
    vector<double> x;
    vector<double> fx;

private:
    solution();
    friend class opt;
};

class opt
{
private:
    opt_func_t f;
    unsigned int d;
    vector<vector<double>> points;
    void *context;

    void sort_by_opt_function();
    vector<double> get_centroid();
    bool reflect(const vector<double> &centroid, double f_1, double f_n, double &f_r, vector<double> &x_r);
    void expand(const vector<double> &centroid, const vector<double> &x_r, double f_r);
    vector<double> outside_contract(const vector<double> &centroid, const vector<double> &x_r);
    bool inside_contract(const vector<double> &centroid, double f_n1);
    void shrink(const vector<double> &x1);
    solution *solution();

public:
    opt(opt_func_t f, void *context, unsigned int d);
    void set_polytope(vector<vector<double>> &points);
    void step();
    void print_points();
    void set_context(void *context);
};

#endif