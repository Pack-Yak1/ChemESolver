#include <stdio.h>
#include <math.h>
#include "../utils/units.h"
using namespace std;

enum log_type
{
    LN,
    LOG
};

class AntoineModel
{
public:
    int A;
    int B;
    int C;
    P_unit p_unit;
    T_unit t_unit;
    log_type log;

    AntoineModel(int a, int b, int c, P_unit p = P_unit::Pa, T_unit t = T_unit::K, log_type l = log_type::LOG)
    {
        A = a;
        B = b;
        C = c;
        p_unit = p;
        t_unit = t;
        log = l;
    }

    // Returns the saturation pressure in units of `p` at temeprature `T`.
    double psat(double T, P_unit p)
    {
        double p = A - B / (C + T);
        p = 
    }
};

int main()
{
    return 0;
}
