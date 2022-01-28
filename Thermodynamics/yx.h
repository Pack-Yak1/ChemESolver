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
    ModifiedRaoultModel(AntoineModel a1, AntoineModel a2, BinaryWilsonModel b,
                        T_unit t_unit, double P, P_unit p_unit);

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

    // Populates `Tx_data` and `Ty_data` with points, where the x value of each
    // point is mole fraction in corresponding liquid/vapor phases, and the y
    // value of each point is temperature.
    // Requires: `num_points` is at least 2, `x_data` and `y_data` are empty,
    // `start` and `end` are between 0 and 1 inclusive and `start` < `end`.
    void generate_Txy_data(int num_points, vector<Point> &Tx_data,
                           vector<Point> &Ty_data, double start = 0,
                           double end = 1);

    // Writes x, y, T values to `o` with `delimiter` between x, y, and T, and
    // `line_break` between each data point.
    // Requires: `num_points` is at least 2, `x_data` and `y_data` are empty,
    // `start` and `end` are between 0 and 1 inclusive and `start` < `end`.
    void write_Txy_data(int num_points, ostream &o = cout, string delim = ",",
                        string line_break = "\n", double start = 0,
                        double end = 1);

    // Helper function which returns the psat of a component modeled by a1 or a2
    // in units of `this->p_unit`, which handles unit conversion of temperature
    // from `this->t_unit`.
    double psat_helper(AntoineModel a, double T);

    // Helper function which returns the tsat of a component modeled by a1 or a2
    // in units of `this->t_unit`, which handles unit conversion of pressire
    // from `this->p_unit`.
    double tsat_helper(AntoineModel a, double P);
};

#endif