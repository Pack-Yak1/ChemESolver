#include <stdlib.h>

#include <iostream>
#include <limits>
#include <vector>

using namespace std;

#ifndef coords
#define coords

class Point
{
public:
    double x;
    double y;

    Point(double x, double y);

    friend ostream &operator<<(std::ostream &os, Point const &p);
};

class Line
{
public:
    double x1;
    double x2;
    double y1;
    double y2;
    double gradient;

    enum error_codes
    {
        INTERSECTION_OF_SAME_LINE = 0,
        INTERSECTION_OF_PARALLEL_LINE,
        LINE_OF_TWO_IDENTICAL_POINTS
    };

    Line(double x1, double y1, double x2, double y2);

    friend ostream &operator<<(std::ostream &os, Line const &l);

    double y_intercept();

    double x_intercept();

    // Returns the point of intersection between the lines:
    // y = `line1[0]` * x + `line1[1]` and
    // y = `line2[0]` * x + `line2[1]`
    static Point intersection(vector<double> &line1, vector<double> &line2);

    static Point intersection(Line l1, Line l2);
};

#endif