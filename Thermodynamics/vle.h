#include "../Thermodynamics/antoine.h"
#include "../Thermodynamics/wilson.h"
#include "../utils/units.h"
#include "../utils/coords.h"
#include <vector>

using namespace std;

class VLE_Interface
{
public:
    virtual void set_Tx(double &x1, double y1, double &T);
    virtual void set_Ty(double x1, double &y1, double &T);
    virtual double psat_helper(AntoineModel a, double T);
    virtual double tsat_helper(AntoineModel a, double P);
    virtual void generate_Txy_data(int num_points, vector<Point> &Tx_data,
                                   vector<Point> &Ty_data, double start = 0,
                                   double end = 1);
    virtual void write_Txy_data(int num_points, ostream &o,
                                string delim, string line_break,
                                double start, double end);
};