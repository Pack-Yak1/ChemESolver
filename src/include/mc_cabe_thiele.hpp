#include "antoine.hpp"
#include "wilson.hpp"
#include "yx.hpp"
#include "units.hpp"
#include <stdbool.h>

#ifndef mc_cabe_thiele
#define mc_cabe_thiele

class MT
{
public:
    double xD;
    double xB;
    double xF;
    bool total_reflux;
    ModifiedRaoultModel m;

    // Returns an instance of a McCabe-Thiele analysis with disillate mole
    // fraction of `xD`, bottoms mole fraction `xB`, feed mole fraction `xF`,
    // which is under total reflux conditions iff `total_reflux` is true.
    MT(ModifiedRaoultModel m, double xD, double xB, double xF, bool total_reflux);

    // Returns a vector with the gradient and y-intercept of the rectifying
    // section operating line as the first and second elements, given the
    // reflux ratio `R`.
    vector<double> rectifying_line(double R);

    // Returns a vector with the gradient and y-intercept of the rectifying
    // section operating line as the first and second elements, given the
    // rectifying section liquid molar flow rate `L` and distillate molar flow
    // rate `D`
    vector<double> rectifying_line(double L, double D);

    // Returns a vector with the gradient and y-intercept of the stripping
    // section operating line as the first and second elements, given the
    // boilup ratio `VB`.
    vector<double> stripping_line(double VB);

    // Returns a vector with the gradient and y-intercept of the stripping
    // section operating line as the first and second elements, given the
    // stripping section vapor molar flow rate `Vbar` and bottoms molar flow
    // rate `B`
    vector<double> stripping_line(double Vbar, double B);

    // Returns the minimum reflux ratio of this instance of McCabe-Thiele
    // analysis assuming saturated liquid feed (i.e. q == 1) and a pinch about
    // the feed.
    double min_reflux();

    // Returns the minimum reflux ratio of this instance of McCabe-Thiele
    // analysis for the value of `q`, assuming a pinch about the feed.
    double min_reflux(double q);

    // Returns a vector of points plotting a pseudo-equilibrium curve from
    // `this->xB` to `this->xD` with size `num_points` for the given
    // `efficiency`.
    // `rect_line` is a vector with the gradient and y-intercept of the
    // rectifying section operating line as the first and second elements.
    // `strip_line` is a vector with the gradient and y-intercept
    // of the stripping section operating line as the first and second elements
    vector<Point> pseudo_equilibrium_curve(unsigned int num_points,
                                           double efficiency,
                                           vector<double> rect_line,
                                           vector<double> strip_line,
                                           int n_workers);

    // Returns a vector of points plotting a pseudo-equilibrium curve from
    // `this->xB` to `this->xD` with size `num_points` for the given
    // `efficiency`. `R` is the reflux ratio. `VB` is the boilup ratio.
    vector<Point> pseudo_equilibrium_curve(unsigned int num_points,
                                           double efficiency,
                                           double R,
                                           double VB,
                                           int n_workers);

    // Returns the number of stages needed for a distillation column under
    // total reflux to separate feed with mole fraction `this->xF` to bottoms
    // with `this->xB` and distillate with `this->xD`. If `verbose` is true,
    // the equilibrium stage mole fractions will be written to `o` with
    // `delim` between mole fractions and `line_break` between each stage.
    int stage_count(bool verbose = false, ostream &o = cout, string delim = ",",
                    string line_break = "\n");

    // Returns the number of stages needed for a distillation column with
    // rectifying and stripping lines described by `rect_line` and `strip_line`
    // to separate feed with mole fraction `this->xF` to bottoms with
    // `this->xB` and distillate with `this->xD`. If `verbose` is true, the
    // equilibrium stage mole fractions will be written to `o` with `delim`
    // between mole fractions and `line_break` between each stage.
    int stage_count(vector<double> rect_line, vector<double> strip_line,
                    bool verbose = false, ostream &o = cout, string delim = ",",
                    string line_break = "\n");
};

#endif