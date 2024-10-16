// https://github.com/NJdevPro/Numerical-Recipes

#pragma once
#include <BSMPT/baryo_calculation/transport_network.h>
#include <BSMPT/utility/matrix_operations.h>
#include <iostream>
#include <memory>

namespace BSMPT
{
class Difeq
{
private:
  // List of z coordinates on the grid
  const VecDoub &z;
  // Equation network we want to solve
  TransportNetwork TN;
  // Jacobian elements (Ainv * (Collision - B))
  Mat3DDoub M;
  // Source term elements (Ainv * Source)
  MatDoub Stilde;

public:
  Difeq(
      // z list
      VecDoub &z_In,
      TransportNetwork tn_in);

  void smatrix(const int k,
               const int k1,
               const int k2,
               const int jsf,
               const int is1,
               const int isf,
               VecInt &indexv,
               MatDoub &s,
               MatDoub &y);
};

} // namespace BSMPT
