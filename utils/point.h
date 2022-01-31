#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits>

using namespace std;

#ifndef point
#define point

class Point
{
public:
    double x;
    double y;

    Point(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    friend ostream &operator<<(std::ostream &os, Point const &p)
    {
        return os << '(' << p.x << ", " << p.y << ')';
    }
};

class Line
{
public:
    double x1;
    double x2;
    double y1;
    double y2;
    double gradient;

    Line(double x1, double y1, double x2, double y2)
    {
        if (x1 == x2 && y1 == y2)
        {
            throw LINE_OF_TWO_IDENTICAL_POINTS;
        }
        this->x1 = x1;
        this->x2 = x2;
        this->y1 = y1;
        this->y2 = y2;
        this->gradient = x1 == x2 ? numeric_limits<double>::infinity() : (y2 - y1) / (x2 - x1);
    }

    friend ostream &operator<<(std::ostream &os, Line const &l)
    {
        return os << "[(" << l.x1 << ", " << l.y1 << "), (" << l.x2 << ", " << l.y2 << ")]";
    }

    double y_intercept()
    {
        return this->y1 - this->x1 * gradient;
    }

    double x_intercept()
    {
        return this->x1 - this->y1 / gradient;
    }

    // Returns the point of intersection between the lines:
    // y = `line1[0]` * x + `line1[1]` and
    // y = `line2[0]` * x + `line2[1]`
    static Point intersection(vector<double> line1, vector<double> line2)
    {
        if (line1[0] == line2[0])
        {
            if (line1[1] == line2[1])
            {
                throw INTERSECTION_OF_SAME_LINE;
            }
            else
            {
                throw INTERSECTION_OF_PARALLEL_LINE;
            }
        }
        double a1 = line1[0];
        double b1 = -1;
        double c1 = line1[1];
        double a2 = line2[0];
        double b2 = -1;
        double c2 = line2[1];
        return Point((b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1),
                     (a2 * c1 - a1 * c2) / (a1 * b2 - a2 * b1));
    }

    static Point intersection(Line l1, Line l2)
    {
        vector<double> line1;
        line1.reserve(2);
        vector<double> line2;
        line2.reserve(2);
        line1.emplace_back(l1.gradient);
        line2.emplace_back(l2.gradient);
        line1.emplace_back(l1.y_intercept());
        line2.emplace_back(l2.y_intercept());
        return intersection(line1, line2);
    }

    enum error_codes
    {
        INTERSECTION_OF_SAME_LINE = 0,
        INTERSECTION_OF_PARALLEL_LINE,
        LINE_OF_TWO_IDENTICAL_POINTS
    };
};

#endif