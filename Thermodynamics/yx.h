#include "antoine.h"
#include "wilson.h"
#include "vector"

#ifndef yx
#define yx

const int NUM_POINTS = 100;

class ModifiedRaoultModel
{
public:
    AntoineModel a1;
    AntoineModel a2;
    BinaryWilsonModel b;
    T_unit t_unit;
    double P;
    P_unit p_unit;

    // An instance of a ModifiedRaoultModel for a binary component system.
    // Relies on Antoine models for each component and a Wilson model for liquid
    // activities.
    ModifiedRaoultModel(AntoineModel a1, AntoineModel a2, BinaryWilsonModel b, T_unit t_unit, double P, P_unit p_unit);

    // Finds a value of mole fraction of component 1 which satisfies the modified
    // Raoult model for temperature `T` in K.
    double find_x1(double T);

    // Finds a value of mole fraction of component 1 which satisfies the modified
    // Raoult model for temperature `T` in units of `t_unit`.
    double find_x1(double T, T_unit t_unit);
};

#endif