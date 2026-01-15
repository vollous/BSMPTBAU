// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/utility/relaxation/difeq.h>

using BSMPT::Delta;

struct Difeq_TransportEquation : Difeq
{
  // List of z coordinates on the grid
  std::vector<double> &z;
  // Number of fermion and number of bosons
  const size_t &nFermions, &nBosons;
  // Tensor with  A^-1 Gamma for each k
  const Mat3DDoub &MTilde;
  // Matrix Stilde = A^-1 * S (source) for each k
  const MatDoub &STilde;
  // Moments used to solve the ODE
  const size_t &moment;

  Difeq_TransportEquation(
      // z list
      std::vector<double> &z_In,
      // Number of fermions and bosons
      const size_t &nFermions_In,
      const size_t &nBosons_In,
      const Mat3DDoub &MTilde_In,
      const MatDoub &STilde_In,
      const size_t &moment_In)
      : z(z_In)
      , nFermions(nFermions_In)
      , nBosons(nBosons_In)
      , MTilde(MTilde_In) // A^-1 * Gamma matrix
      , STilde(STilde_In) // A^-1 * Sources vector
      , moment(moment_In)
  {
  }

  void smatrix(const int k,
               const int k1,
               const int k2,
               const int jsf,
               VecInt &indexv,
               MatDoub &s,
               MatDoub &y)
  {
    const size_t nP = nFermions + nBosons;
    double temp;
    s.zero(); // Set matrix s = 0
    if (k == k1)
    {
      // Boundary conditions mu = 0 on first boundary
      for (size_t particle = 0; particle < nP; particle++)
      {
        // Sn at the first boundary
        s[nP + particle][moment * nP + indexv[moment * particle]] = 1.0;
        // B0
        s[nP + particle][jsf] = y[indexv[moment * particle]][0];
      }
    }
    else if (k > k2 - 1)
    {
      // Boundary conditions mu = 0 on second boundary
      for (size_t particle = 0; particle < nP; particle++)
      {
        // Sn at the last boundary
        s[particle][moment * nP + indexv[moment * particle]] = 1.0;
        // C0
        s[particle][jsf] = y[indexv[moment * particle]][z.size() - 1];
      }
    }
    else
    {
      for (size_t j = 0; j < moment * (nFermions + nBosons); j++)
      {
        for (size_t n = 0; n < moment * (nFermions + nBosons); n++)
        {
          // s matrix for the middle point
          // S_{j,n}
          // M[k] and S[k] are evaluated between k and k-1
          s[j][indexv[n]] =
              -Delta(j, n) - 0.5 * (z[k] - z[k - 1]) * (MTilde[k][j][n]);
          // S_{j,N + n}
          s[j][2 * nP + indexv[n]] =
              Delta(j, n) - 0.5 * (z[k] - z[k - 1]) * (MTilde[k][j][n]);
        }
        //  Equations for E(k,k-1)
        temp = STilde[k][j];
        for (size_t i = 0; i < moment * (nFermions + nBosons); i++)
        {
          temp += MTilde[k][j][i] * (y[i][k] + y[i][k - 1]) / 2.;
        }
        s[j][jsf] = y[j][k] - y[j][k - 1] - (z[k] - z[k - 1]) * temp;
      }
    }
  }
};
