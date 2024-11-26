// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/utility/data_structures.h>
#include <BSMPT/utility/utility.h>

using BSMPT::Delta;

struct Difeq
{
  // List of z coordinates on the grid
  VecDoub &z;
  // Number of fermion and number of bosons
  const int &nFermions, &nBosons;
  // Tensor with  A^-1 Gamma for each k
  const Mat3DDoub &MTilde;
  // Matrix Stilde = A^-1 * S (source) for each k
  const MatDoub &STilde;

  Difeq(
      // z list
      VecDoub &z_In,
      // Number of fermions and bosons
      const int &nFermions_In,
      const int &nBosons_In,
      const Mat3DDoub &MTilde_In,
      const MatDoub &STilde_In)
      : z(z_In)
      , nFermions(nFermions_In)
      , nBosons(nBosons_In)
      , MTilde(MTilde_In) // A^-1 * Gamma matrix
      , STilde(STilde_In) // A^-1 * Sources vector
  {
  }

  void smatrix(const int k,
               const int k1,
               const int k2,
               const int jsf,
               const int is1,
               const int isf,
               VecInt &indexv,
               MatDoub &s,
               MatDoub &y)
  {
    (void)is1;
    (void)isf;
    const int nP = nFermions + nBosons;
    double temp;
    s.zero(); // Set matrix s = 0
    if (k == k1)
    {
      // Boundary conditions mu = 0 on first boundary
      for (int particle = 0; particle < nP; particle++)
      {
        // Sn at the first boundary
        s[nP + particle][2 * nP + indexv[2 * particle]] = 1.0;
        // B0
        s[nP + particle][jsf] = y[indexv[2 * particle]][0];
      }
    }
    else if (k > k2 - 1)
    {
      // Boundary conditions mu = 0 on second boundary
      for (int particle = 0; particle < nP; particle++)
      {
        // Sn at the last boundary
        s[particle][2 * nP + indexv[2 * particle]] = 1.0;
        // C0
        s[particle][jsf] = y[indexv[2 * particle]][z.size() - 1];
      }
    }
    else
    {
      for (int j = 0; j < 2 * (nFermions + nBosons); j++)
      {
        for (int n = 0; n < 2 * (nFermions + nBosons); n++)
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
        for (int i = 0; i < 2 * (nFermions + nBosons); i++)
        {
          temp += MTilde[k][j][i] * (y[i][k] + y[i][k - 1]) / 2.;
        }
        s[j][jsf] = y[j][k] - y[j][k - 1] - (z[k] - z[k - 1]) * temp;
      }
    }
  }
};
