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
    vector<vector<double>> points;

    void sort_by_opt_function();
    vector<double> get_centroid();
    vector<double> reflect(const vector<double> &centroid, double f_1, double f_n, double &f_r);
    void expand(const vector<double> &centroid, const vector<double> &x_r, double f_r);
    vector<double> outside_contract(const vector<double> &centroid, const vector<double> &x_r);
    void inside_contract(const vector<double> &centroid, const vector<double> &x_r, double f_n1);
    void shrink(const vector<double> &x1);
    
public:
    opt(opt_func_t f, unsigned int d);
    void set_polytope(vector<vector<double>> &points);
    void step();
};

#endif