#include "mc_cabe_thiele.hpp"

#include <vector>

#include "antoine.hpp"
#include "wilson.hpp"
#include "yx.hpp"
#include "coords.hpp"
#include "units.hpp"
#include "opt.hpp"

using namespace std;

/**
 * @brief Class for passing a linear constraint and model data to opt solver
 */
class opt_context
{
public:
    ModifiedRaoultModel *m;
    double grad;
    double y_intersect;

    opt_context(ModifiedRaoultModel *m, double grad, double y_intersect)
    {
        this->m = m;
        this->grad = grad;
        // 0 for x1, 1 for y1, 2 for T.
        this->y_intersect = y_intersect;
    }
};

MT::MT(ModifiedRaoultModel m, double xD, double xB, double xF,
       bool total_reflux)
{
    this->m = m;
    this->xD = xD;
    this->xB = xB;
    this->xF = xF;
    this->total_reflux = total_reflux;
}

vector<double> MT::rectifying_line(double R)
{
    vector<double> output;
    output.reserve(2);
    output.emplace_back(R / (R + 1));
    output.emplace_back(xD / (R + 1));
    // cout << output[0] << "," << output[1] << "\n";
    return output;
}

vector<double> MT::rectifying_line(double L, double D)
{
    double R = L / D;
    return rectifying_line(R);
}

vector<double> MT::stripping_line(double VB)
{
    vector<double> output;
    output.reserve(2);
    output.emplace_back((VB + 1) / VB);
    output.emplace_back(-xB / VB);
    // cout << output[0] << "," << output[1] << "\n";
    return output;
}

vector<double> MT::stripping_line(double Vbar, double B)
{
    double VB = Vbar / B;
    return stripping_line(VB);
}

// Assumes saturated liquid feed
double MT::min_reflux()
{
    double yF = xF;
    double T;
    m.set_Ty(xF, yF, T);
    // cout << xF << ',' << yF << ',' << T << '\n';
    double gradient = Line(xF, yF, xD, xD).gradient;
    return -gradient / (gradient - 1);
}

double min_reflux_objective(const vector<double> &x, void *context)
{
    opt_context *ctx = (opt_context *)context;
    ModifiedRaoultModel *m = ctx->m;
    double grad = ctx->grad;
    double its = ctx->y_intersect;

    double x1 = x[0];
    double y1 = x[1];
    double T = x[2];

    double gamma1 = m->b.gamma1(x1, T, m->t_unit);
    double gamma2 = m->b.gamma2(1 - x1, T, m->t_unit);
    double psat1 = m->psat_helper(m->a1, T);
    double psat2 = m->psat_helper(m->a2, T);
    double fugacity1 = (gamma1 * x1 * psat1);
    double fugacity2 = (gamma2 * (1 - x1) * psat2);

    double constraint_loss = pow(x1 * grad + its - y1, 2.);
    double y1_loss = pow(fugacity1 / m->P - y1, 2.);
    double y2_loss = pow(fugacity2 / m->P - (1 - y1), 2.);
    return constraint_loss + y1_loss + y2_loss;
}

// Finds minimum reflux ratio for a specified q-value
double MT::min_reflux(double q)
{
    if (q == 1)
    {
        return min_reflux();
    }

    // Formula for q-line from CHEME 3320
    opt_context ctx(&m, -q / (1. - q), -xF / (q - 1.));
    opt solver(min_reflux_objective, &ctx, 3);
    double bp1 = m.tsat_helper(m.a1, m.P);
    double bp2 = m.tsat_helper(m.a2, m.P);
    vector<double> lb{0, 0, 0};
    vector<double> ub{1, 1, 10 * max(bp1, bp2)};

    solution *soln = solver.solve(lb, ub);
    double x_pinch = soln->x[0];
    double y_pinch = soln->x[1];
    double gradient = Line(x_pinch, y_pinch, xD, xD).gradient;
    delete soln;
    return -gradient / (gradient - 1);
}

vector<Point> MT::pseudo_equilibrium_curve(unsigned int num_points,
                                           double efficiency,
                                           vector<double> rect_line,
                                           vector<double> strip_line,
                                           int n_workers)
{
    Point intersection = Line::intersection(rect_line, strip_line);
    vector<Point> Tx_data;
    vector<Point> Ty_data;
    m.generate_Txy_data(num_points, Tx_data, Ty_data, xB, xD, n_workers);
    vector<Point> output;
    output.reserve(num_points);
    for (unsigned int i = 0; i < num_points; i++)
    {
        double x = Tx_data[i].x;
        double y = Ty_data[i].x;
        double y_op;
        if (total_reflux)
        {
            y_op = x;
        }
        else
        {
            if (x < intersection.x)
            {
                y_op = strip_line[0] * x + strip_line[1];
            }
            else
            {
                y_op = rect_line[0] * x + rect_line[1];
            }
        }
        double y_out = (y - y_op) * efficiency + y_op;
        output.emplace_back(x, y_out);
    }
    return output;
}

vector<Point> MT::pseudo_equilibrium_curve(unsigned int num_points,
                                           double efficiency, double R,
                                           double VB, int n_workers)
{
    return pseudo_equilibrium_curve(num_points, efficiency, rectifying_line(R),
                                    stripping_line(VB), n_workers);
}

// For total reflux
int MT::stage_count(bool verbose, ostream &o, string delim, string line_break)
{
    double x = xD;
    double y = xD;
    double T;
    int stage_count = 0;
    o << "Tray Number" << delim << "x" << delim << "y" << line_break;
    while (x > xB)
    {
        m.set_Tx(x, y, T);
        stage_count++;
        if (verbose)
        {
            o << stage_count << delim << x << delim << y << line_break;
        }
        y = x;
    }
    return stage_count;
}

int MT::stage_count(vector<double> rect_line, vector<double> strip_line,
                    bool verbose, ostream &o, string delim, string line_break)
{
    double x = xD;
    double y = xD;
    double T;
    int stage_count = 0;
    o << "Tray Number" << delim << "x" << delim << "y" << line_break;
    Point intersection = Line::intersection(rect_line, strip_line);
    while (x > xB)
    {
        m.set_Tx(x, y, T);
        // printf("stage count: %d, x: %f, y: %f, T: %f\n", stage_count, x, y,
        // T);
        stage_count++;
        if (verbose)
        {
            o << stage_count << delim << x << delim << y << line_break;
        }
        vector<double> *op_line = x > intersection.x ? &rect_line : &strip_line;
        y = (*op_line)[0] * x + (*op_line)[1];
        if (verbose)
        {
            o << stage_count << delim << x << delim << y << line_break;
        }
    }
    return stage_count;
}