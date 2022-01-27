#include "../Thermodynamics/antoine.h"
#include "../Thermodynamics/wilson.h"
#include "../Thermodynamics/yx.h"
#include "../utils/units.h"
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

    MT(ModifiedRaoultModel m, double xD, double xB, double xF, bool total_reflux);

    vector<double> rectifying_line(double L, double D);
    vector<double> stripping_line(double L, double D);
    double min_reflux();
    double min_reflux(double q);
};

#endif