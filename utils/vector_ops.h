#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

#define DIMENSION_ERROR 1
#define DEGENERATE_POINTS 2

vector<double> sum(vector<double> v1, vector<double> v2);
vector<double> sum(const vector<vector<double>> &vecs);
vector<double> linspace(double start, double end, double points);
vector<vector<double>> zip(vector<vector<double>> linspaces);
