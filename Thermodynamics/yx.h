#include "antoine.h"
#include "wilson.h"
#include "vector"
#include "../utils/point.h"
#include "../utils/units.h"
#include "gnuplot-iostream.h"
using namespace std;

#ifndef yx
#define yx

const double REL_XTOL = 1e-10;

double f(const std::vector<double> &x, std::vector<double> &grad, void *f_data);

class ModifiedRaoultModel
{
public:
    AntoineModel a1;
    AntoineModel a2;
    BinaryWilsonModel b;
    T_unit t_unit;
    double P;
    P_unit p_unit;

    ModifiedRaoultModel();

    // An instance of a ModifiedRaoultModel for a binary component system.
    // Relies on Antoine models for each component and a Wilson model for liquid
    // activities.
    ModifiedRaoultModel(AntoineModel a1, AntoineModel a2, BinaryWilsonModel b, T_unit t_unit, double P, P_unit p_unit);

    // Finds a value of mole fraction of component 1 which satisfies the modified
    // Raoult model for temperature `T` in K.
    double find_x1(double T);

    // Finds a value of mole fraction of component 1 which satisfies the modified
    // Raoult model for temperature `T` in units of `t_unit`.
    double find_x1(double T, T_unit t_unit);

    // Finds a value of `T` in units of K that is consistent with this model's
    // pressure and the value of `x1`.
    double find_T(double x1);

    // Sets the values of `T` in units of K and `y` to values consistent with
    // this model's pressure and the value of `x1`
    void find_Ty(double x1, double &y1, double &T);

    void generate_Txy_data(int num_points, vector<Point> &Tx_data, vector<Point> &Ty_data);

    void write_Txy_data(int num_points, ostream &o = cout, string delim = ",", string line_break = "\n");

    double psat_helper(AntoineModel a, double T);

    double tsat_helper(AntoineModel a, double P);
};

#endif