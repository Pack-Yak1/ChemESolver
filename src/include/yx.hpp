#include <mutex>

#include "coords.hpp"
#include "units.hpp"
#include "antoine.hpp"
#include "vector"
#include "wilson.hpp"

using namespace std;

#ifndef yx
#define yx

int const x_var = 0;
int const y_var = 1;
int const T_var = 2;

const double REL_XTOL = 1e-10;

double p_error(const vector<double> &x, void *context);
double x_error(const std::vector<double> &x, void *f_data);

class ModifiedRaoultModel
{
private:
    void generate_Txy_single_thread(double start, double step_size,
                                    double num_steps, int pos,
                                    vector<Point> &Tx_data,
                                    vector<Point> &Ty_data, mutex &mtx);

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

    /**
     * @brief Solve the system for liquid and vapor phase mole fractions of
     * component 1 and temperature, given any of the 3. Denoted as (x1, y1, T),
     * respectively.
     *
     * @param val The known value among x1, y1, T.
     * @param dim Must be x_var if x1, y_var if y1, T_var if T.
     * @return Size 3 vector containing [x1, y1, T].
     */
    vector<double> solve_from_constraint(double val, int dim);

    /**
     * @brief Solve the system for liquid and vapor phase mole fractions of
     * component 1 and temperature in a specified unit, given any of the 3.
     * Denoted as (x1, y1, T), respectively.
     *
     * @param val The known value among x1, y1, T.
     * @param dim Must be x_var if x1, y_var if y1, T_var if T.
     * @param t_unit The unit which T is in. Unused if `dim` is not T_var
     * @return Size 3 vector containing [x1, y1, T].
     */
    vector<double> solve_from_constraint(double val, int dim, T_unit t_unit);

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