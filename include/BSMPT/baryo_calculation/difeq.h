// https://github.com/NJdevPro/Numerical-Recipes

#pragma once
#include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/utility/nr3.h>
#include <BSMPT/utility/utility.h>

namespace BSMPT
{
class Difeq
{
private:
  // Information about vw and Tc
  std::shared_ptr<Kinfo> Ki;
  // List of z coordinates on the grid
  const VecDoub &z;
  // Number of fermion and number of bosons
  const Int &nFermions, &nBosons;
  // Matrix Stilde = A^-1 * S (source) for each k
  MatDoub STilde;
  // Tensor with  A^-1 Gamma for each k
  Mat3DDoub M;

public:
  Difeq(
      std::shared_ptr<Kinfo> K_in,
      // z list
      VecDoub &z_In,
      // Number of fermions and bosons
      const Int &nFermions_In,
      const Int &nBosons_In,
      // Fermions
      const MatDoub &fm2prime,
      const MatDoub &fS,
      // Bosons
      const MatDoub &bm2prime,
      // Collision matrix
      const Mat3DDoub &CollisionMatrix); // A^-1 * Sources vector

  void smatrix(const Int k,
               const Int k1,
               const Int k2,
               const Int jsf,
               const Int is1,
               const Int isf,
               VecInt_I &indexv,
               MatDoub_O &s,
               MatDoub_I &y)
  {
    (void)is1;
    (void)isf;
    Doub temp;
    if (k == k1)
    {
      // Boundary conditions mu = 0 on first boundary
      for (int particle = 0; particle < nFermions + nBosons; particle++)
      {
        // Sn at the first boundary
        s[nFermions + nBosons + particle][indexv[2 * particle]] = 1.0;
        // B0
        s[nFermions + nBosons + particle][jsf] = y[indexv[2 * particle]][0];
      }
    }
    else if (k > k2 - 1)
    {
      // Boundary conditions mu = 0 on second boundary
      for (int particle = 0; particle < nFermions + nBosons; particle++)
      {
        // Sn at the last boundary
        s[particle][indexv[2 * particle]] = 1.0;
        // C0
        s[particle][jsf] = y[indexv[2 * particle]][z.size() - 1];
      }
    }
    else
    {
      for (int j = 0; j < 2 * (nFermions + nBosons); j++)
        for (int n = 0; n < 2 * (nFermions + nBosons); n++)
        {
          // s matrix for the middle point
          // S_{j,n}
          s[j][indexv[n]] =
              -Delta(j, n) - 0.5 * (z[k] - z[k - 1]) * (M[k - 1][j][n]);
          // S_{j,N + n}
          s[j][2 * (nFermions * nBosons) + indexv[n]] =
              Delta(j, n) - 0.5 * (z[k] - z[k - 1]) * (M[k][j][n]);
          // Equations for E(k,k-1)
          temp = STilde[k][j] + STilde[k - 1][j];
          for (int i = 0; i <= 2 * (nFermions + nBosons); i++)
            temp += M[k][j][i] * y[i][k] + M[k - 1][j][i] * y[i][k - 1];
          s[j][jsf] = y[j][k] - y[j][k - 1] - 0.5 * (z[k] - z[k - 1]) * temp;
        }
    }
  }
};

} // namespace BSMPT
