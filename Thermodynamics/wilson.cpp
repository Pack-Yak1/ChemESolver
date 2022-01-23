#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../utils/units.h"
#include <stdbool.h>
using namespace std;

const double IDEAL_GAS_CONST = 8.314;

class BinaryWilsonModel
{
public:
  double v1;
  double v2;
  double delta_l12;
  double delta_l21;
  double L12;
  double L21;
  bool hasTemp;

  // Initializes an instance of BinaryWilsonModel which does not have a
  // specified temperature
  BinaryWilsonModel(double v1, double v2, double delta_l12, double delta_l21)
  {
    this->v1 = v1;
    this->v2 = v2;
    this->delta_l12 = delta_l12;
    this->delta_l21 = delta_l21;
    hasTemp = false;
  }

  // Initializes an instance of BinaryWilsonModel with a specified temperature
  // `T` in units of `t_unit`
  BinaryWilsonModel(double v1, double v2, double delta_l12, double delta_l21, double T, T_unit t_unit)
  {
    this->v1 = v1;
    this->v2 = v2;
    this->delta_l12 = delta_l12;
    this->delta_l21 = delta_l21;
    hasTemp = true;
    double t = convert_T(T, t_unit, T_unit::K);
    L12 = v2 / v1 * exp(-delta_l12 / IDEAL_GAS_CONST / t);
    L21 = v1 / v2 * exp(-delta_l21 / IDEAL_GAS_CONST / t);
  }

  void setT(double T, T_unit t_unit)
  {
    double t = convert_T(T, t_unit, T_unit::K);
    L12 = v2 / v1 * exp(-delta_l12 / IDEAL_GAS_CONST / t);
    L21 = v1 / v2 * exp(-delta_l21 / IDEAL_GAS_CONST / t);
    hasTemp = true;
  }

  // Returns the activity coefficient of component 1 for a BinaryWilsonModel
  // which already has a specified temperature
  double gamma1(double x1)
  {
    if (!hasTemp)
    {
      throw "Called gamma1 from BinaryWilsonModel with no specified temperature.\n";
    }
    double x2 = 1 - x1;
    return exp(
        -log(x1 + x2 * L12) + x2 * (L12 / (x1 + x2 * L12) - L21 / (x1 * L21 + x2)));
  }

  // Destructively modifies this BinaryWilsonModel and returns the acitivity
  // coefficient of component 1 at temperature `T` in units of `t_unit`.
  double gamma1(double x1, double T, T_unit t_unit)
  {
    setT(T, t_unit);
    return gamma1(x1);
  }

  double gamma2(double x2)
  {
    if (!hasTemp)
    {
      throw "Called gamma2 from BinaryWilsonModel with no specified temperature.\n";
    }
    double x1 = 1 - x2;
    return exp(
        -log(x2 + x1 * L21) - x1 * (L12 / (x1 + x2 * L12) - L21 / (x1 * L21 + x2)));
  }

  // Destructively modifies this BinaryWilsonModel and returns the acitivity
  // coefficient of component 2 at temperature `T` in units of `t_unit`.
  double gamma2(double x2, double T, T_unit t_unit)
  {
    setT(T, t_unit);
    return gamma2(x2);
  }
};

int main()
{
  return 0;
}
