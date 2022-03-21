#include "../Thermodynamics/antoine.h"
#include "../Thermodynamics/wilson.h"
#include "../Thermodynamics/yx.h"
#include "../utils/units.h"
#include "../utils/coords.h"
#include "../Distillation/mc_cabe_thiele.h"
#include <nlopt.hpp>
#include <vector>

using namespace std;

int main()
{
    // TODO: add gui
    // TODO: use references to pass args
    // TODO: consider making lower level solver methods private
    AntoineModel a1 = antoine::ANTOINE_METHANOL;
    AntoineModel a2 = antoine::ANTOINE_WATER;

    BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398);
    ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);
    MT mt(m, 0.98, 0.02, 0.2, false);

    vector<Point> peq = mt.pseudo_equilibrium_curve(101, 0.1, 1.86805487, 0.829084445);
    cout << "Pseudoequilibrium curve generated, compare against lab data:\n";

    for (auto i = peq.begin(); i != peq.end(); i++)
    {
        cout << *i << "\n";
    }
    // vector<double> v;
    // v.reserve(2);
    // v.emplace_back(1);
    // v.emplace_back(0);

    ofstream output_file;
    output_file.open("seps.csv");
    // mt.m.write_Txy_data(1001, output_file);
    mt.m.write_Txy_data(101, output_file);

    // for (int i = 0; i < peq.size(); i++)
    // {
    //     cout << peq[i].x << "," << peq[i].y << '\n';
    // }

    // double x = .8;
    // double y = .91;
    // double T = 358;
    // mt.m.set_xy(x, y, T);
    // printf("x: %f, y: %f, T: %f\n", x, y, T);
    // cout << mt.m.solve_from_y1(.78) << '\n';
    // cout << mt.m.find_T(.78, false) << '\n';
    // cout << mt.m.find_T(0.4859, true) << '\n';
    // cout << mt.stage_count(mt.rectifying_line(1.04027 * 1.3), mt.stripping_line(0.542850835), true, output_file) << '\n';
    // mt.stage_count(true, output_file, ",", "\n");
    // Line l1 = Line(0, 0, 1, 1);
    // Line l2 = Line(1, 1.5, 2.5, 0);
    // cout << Line::intersection(l1, l2);

    // vector<double> rect = mt.rectifying_line(mt.min_reflux() * 1);
    // cout << rect[0] << ',' << rect[1];
    return 0;
}