#include <iostream>
#include <stdlib.h>

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

    friend std::ostream &operator<<(std::ostream &os, Point const &p)
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

    Line(double x1, double y1, double x2, double y2)
    {
        this->x1 = x1;
        this->x2 = x2;
        this->y1 = y1;
        this->y2 = y2;
    }

    friend std::ostream &operator<<(std::ostream &os, Line const &l)
    {
        return os << "[(" << l.x1 << ", " << l.y1 << "), (" << l.x2 << ", " << l.y2 << ")]";
    }

    double gradient()
    {
        return (y2 - y1) / (x2 - x1);
    }
};

#endif