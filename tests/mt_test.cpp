#include <nlopt.hpp>
#include <vector>

#include "../Distillation/mc_cabe_thiele.h"
#include "../Thermodynamics/antoine.h"
#include "../Thermodynamics/wilson.h"
#include "../Thermodynamics/yx.h"
#include "../utils/coords.h"
#include "../utils/units.h"

using namespace std;

const int NUM_POINTS = 100001;
const int MAX_THREADS = 100;

int main()
{
    AntoineModel a1 = antoine::ANTOINE_METHANOL;
    AntoineModel a2 = antoine::ANTOINE_WATER;

    BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398);
    ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);
    MT mt(m, 0.98, 0.02, 0.2, false);

    vector<Point> peq =
        mt.pseudo_equilibrium_curve(1001, 0.1, 1.86805487, 0.829084445, 10);
    cout << "Pseudoequilibrium curve generated, compare against lab data:\n";

    for (auto i = peq.begin(); i != peq.end(); i++)
    {
        cout << *i << "\n";
    }

    // ofstream output_file;
    // output_file.open("seps.csv");
    // mt.m.write_Txy_data(10001, output_file);
    cout << "Test time taken under different number of workers\n";
    for (int i = 1; i < MAX_THREADS; i++)
    {
        vector<Point> Tx_data;
        vector<Point> Ty_data;
        auto start = chrono::high_resolution_clock::now();
        mt.m.generate_Txy_data(NUM_POINTS, Tx_data, Ty_data, 0, 1, i);
        auto stop = chrono::high_resolution_clock::now();
        double duration = chrono::duration<double>(stop - start).count();
        cout << "Time taken to calculate " << NUM_POINTS << " points with " << i
             << " worker threads was " << duration << "s.\n";
    }
    return 0;
}