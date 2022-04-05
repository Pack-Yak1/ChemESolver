#include <vector>
#include <thread>

#include "mc_cabe_thiele.hpp"
#include "antoine.hpp"
#include "wilson.hpp"
#include "yx.hpp"
#include "coords.hpp"
#include "units.hpp"

using namespace std;

const int NUM_POINTS = 10;
const unsigned int MAX_THREADS = thread::hardware_concurrency();

int main()
{
    AntoineModel a1 = antoine::ANTOINE_METHANOL;
    AntoineModel a2 = antoine::ANTOINE_WATER;

    BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398);
    ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);
    MT mt(m, 0.98, 0.02, 0.2, false);

    vector<Point> peq = mt.pseudo_equilibrium_curve(
        NUM_POINTS, 0.1, 1.86805487, 0.829084445, MAX_THREADS);
    cout << "Pseudoequilibrium curve generated, compare against lab data:\n";

    for (auto i = peq.begin(); i != peq.end(); i++)
    {
        cout << *i << "\n";
    }

    cout << "Testing min reflux, compare against continuous distillation data\n";
    cout << "For saturated liquid, R_min is: " << mt.min_reflux() << "\n";
    for (double q = -1.; q < 0.; q += 0.1)
    {
        cout << "For q = " << q << ", R_min is: " << mt.min_reflux(q) << "\n";
    }
    // cout << "For saturated vapor, R_min is: " << mt.min_reflux(0) << "\n";
    // cout << "For q = 0.5, R_min is: " << mt.min_reflux(0.5) << "\n";

    // ofstream output_file;
    // output_file.open("seps.csv");
    // mt.m.write_Txy_data(10001, output_file);
}