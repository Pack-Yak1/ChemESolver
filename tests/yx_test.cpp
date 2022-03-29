#include <vector>
#include <thread>

#include "../Thermodynamics/antoine.h"
#include "../Thermodynamics/wilson.h"
#include "../Thermodynamics/yx.h"
#include "../utils/coords.h"
#include "../utils/units.h"

using namespace std;

const int NUM_POINTS = 1000000;
const unsigned int NUM_CORES = thread::hardware_concurrency();

void test_multithreaded_Ty_solver(ModifiedRaoultModel m)
{
    for (unsigned int i = 1; i <= NUM_CORES; i++)
    {
        // Tests solving system from x1
        ofstream output_file;
        string filename("test_yx_from_x1_" + to_string(i) + "_threads" + ".csv");
        output_file.open(filename);
        auto start = chrono::high_resolution_clock::now();
        m.write_Txy_data(NUM_POINTS, output_file, ",", "\n", 0, 1, i);
        auto stop = chrono::high_resolution_clock::now();
        double duration = chrono::duration<double>(stop - start).count();
        output_file.close();
        cout << "Time taken to calculate " << NUM_POINTS << " points with " << i
             << " worker threads was " << duration << "s.\n";
        cout << "Data sent to " << filename << ". Verify against known data.\n\n";
    }
}

void set_xy_test(ModifiedRaoultModel m, AntoineModel a1, AntoineModel a2)
{
    double x1, y1;
    ofstream output_file;
    output_file.open("set_xy_test.txt");
    for (double temp = m.tsat_helper(a1, m.P); temp < m.tsat_helper(a2, m.P); temp += 0.1)
    {
        m.set_xy(x1, y1, temp);
        output_file << "At T = " << temp
                    << "K, the equilibrium composition is ("
                    << x1 << ", " << y1 << ")\n";
    }
    output_file.close();
}

void set_Tx_test(ModifiedRaoultModel m, AntoineModel a1, AntoineModel a2)
{
    double x1, temp;
    ofstream output_file;
    output_file.open("set_Ty_test.txt");
    for (double y1 = 0.; y1 <= 1; y1 += 0.001)
    {
        m.set_Tx(x1, y1, temp);
        output_file << "At y = " << y1 << ", x = " << x1 << ", T = " << temp
                    << "K\n";
    }
    output_file.close();
}

int main()
{
    AntoineModel a1 = antoine::ANTOINE_METHANOL;
    AntoineModel a2 = antoine::ANTOINE_WATER;

    BinaryWilsonModel b(0.00004073, 0.00001807, 347.4525, 2179.8398);
    ModifiedRaoultModel m(a1, a2, b, T_unit::K, 99250, P_unit::Pa);

    // Test multithreading and solving system given liq mole frac constraint
    cout << "Test time taken under different number of workers\n";
    test_multithreaded_Ty_solver(m);

    // Tests solving system from temperature
    cout << "Test `set_xy`\n";
    // set_xy_test(m, a1, a2);

    // Tests solving system from vap mole frac constraint
    cout << "Test `set_Tx`\n";
    // set_Tx_test(m, a1, a2);

    return 0;
}