#include "../utils/coords.h"
#include "../utils/units.h"
#include "antoine.h"
#include "gnuplot-iostream.h"
#include "vector"
#include "wilson.h"
using namespace std;

#ifndef yx
#define yx

const double REL_XTOL = 1e-10;

double p_error(const vector<double> &x, void *context);
double p_error(const std::vector<double> &x, std::vector<double> &grad,
               void *f_data);
double x_error(const std::vector<double> &x, void *f_data);

class ModifiedRaoultModel
{
private:
    void generate_Txy_single_thread(double start, double step_size,
                                    double num_steps, int pos,
                                    vector<Point> &Tx_data,
                                    vector<Point> &Ty_data);

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

    /**
     * @brief Copy constructor for ModifiedRaoultModel
     */
    ModifiedRaoultModel(const ModifiedRaoultModel &m);

    // Finds mole fractions in liquid and vapor phase consistent with `T` in
    // K. Returns a vector with first element = x1, second = y1, third = T.
    vector<double> solve_from_T(double T);

    // Finds mole fractions in liquid and vapor phase consistent with `T` in
    // units of `t_unit`.
    // Returns a vector with first element = x1, second = y1, third = T.
    vector<double> solve_from_T(double T, T_unit t_unit);

    // Finds mole fraction in liquid phase and temperature consistent with `y1`.
    // Returns a vector with first element = x1, second = y1, third = T.
    vector<double> solve_from_y1(double y1);

    // Finds mole fraction in vapor phase and temperature consistent with `x1`.
    // Returns a vector with first element = x1, second = y1, third = T.
    vector<double> solve_from_x1(double x1);

    // Sets the values of `T` in units of K and `y1` to values consistent with
    // this model's pressure and the value of `x1`
    void set_Ty(double x1, double &y1, double &T);

    // Sets the values of `T` in units of K and `x1` to values consistent with
    // this model's pressure and the value of `y1`
    void set_Tx(double &x1, double y1, double &T);

    // Sets the values of `x1` and `y1` to values consistent with
    // this model's pressure and the value of `T` in `K`.
    void set_xy(double &x1, double &y1, double T);

    // Populates `Tx_data` and `Ty_data` with points, where the x value of each
    // point is mole fraction in corresponding liquid/vapor phases, and the y
    // value of each point is temperature.
    // Requires: `num_points` is at least 2, `x_data` and `y_data` are empty,
    // `start` and `end` are between 0 and 1 inclusive and `start` < `end`.
    void generate_Txy_data(int num_points, vector<Point> &Tx_data,
                           vector<Point> &Ty_data, double start = 0,
                           double end = 1, int n_workers = 4);

    // Writes x, y, T values to `o` with `delimiter` between x, y, and T, and
    // `line_break` between each data point.
    // Requires: `num_points` is at least 2, `x_data` and `y_data` are empty,
    // `start` and `end` are between 0 and 1 inclusive and `start` < `end`.
    void write_Txy_data(int num_points, ostream &o = cout, string delim = ",",
                        string line_break = "\n", double start = 0,
                        double end = 1, int n_workers = 4);

    // Helper function which returns the psat of a component modeled by a1 or a2
    // in units of `this->p_unit`, which handles unit conversion of temperature
    // from `this->t_unit`.
    double psat_helper(AntoineModel a, double T);

    // Helper function which returns the tsat of a component modeled by a1 or a2
    // in units of `this->t_unit`, which handles unit conversion of pressure
    // from `this->p_unit`.
    double tsat_helper(AntoineModel a, double P);
};

#endif