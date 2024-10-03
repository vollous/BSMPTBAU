// https://github.com/NJdevPro/Numerical-Recipes

#pragma once
#include <BSMPT/utility/nr3.h>
#include <BSMPT/utility/utility.h>

using BSMPT::Delta;

struct Difeq
{
  // vw
  const Doub &vw;
  Doub gamma;
  // List of z coordinates on the grid
  const VecDoub &z;
  // Number of fermion and number of bosons
  const Int &nFermions, &nBosons;
  // Thermal averages functions for fermions. Matrices k-steps x nFermion
  const MatDoub &fD0, &fD1, &fD2, &fR, &fRbar, &fQ1, &fQ2, &fm2prime, &fS;
  // Thermal averages functions for bosons. Matrices k-steps x nBoson
  const MatDoub &bD1, &bD2, &bR, &bRbar, &bQ1, &bQ2, &bm2prime;
  // Collision matrix.
  // Rank-3 "tensor" k-steps x 2 * (nFermions + nBosons) x 2 * (nFermions +
  // nBosons)
  const Mat3DDoub &CollisionMatrix;
  // Matrix Stilde = A^-1 * S (source) for each k
  MatDoub STilde;
  // Tensor with  A^-1 Gamma for each k
  Mat3DDoub M;

  Difeq(
      // Wall velocity
      const Doub &vw_In,
      // z list
      VecDoub &z_In,
      // Number of fermions and bosons
      const Int &nFermions_In,
      const Int &nBosons_In,
      // Fermions
      const MatDoub &fD0_In,
      const MatDoub &fD1_In,
      const MatDoub &fD2_In,
      const MatDoub &fR_In,
      const MatDoub &fRbar_In,
      const MatDoub &fQ1_In,
      const MatDoub &fQ2_In,
      const MatDoub &fm2prime_In,
      const MatDoub &fS_In,
      // Bosons
      const MatDoub &bD1_In,
      const MatDoub &bD2_In,
      const MatDoub &bR_In,
      const MatDoub &bRbar_In,
      const MatDoub &bQ1_In,
      const MatDoub &bQ2_In,
      const MatDoub &bm2prime_In,
      // Collision matrix
      const Mat3DDoub &CollisionMatrix_In)
      : vw(vw_In)
      , z(z_In)
      , nFermions(nFermions_In)
      , nBosons(nBosons_In)
      , fD0(fD0_In)
      , fD1(fD1_In)
      , fD2(fD2_In)
      , fR(fR_In)
      , fRbar(fRbar_In)
      , fQ1(fQ1_In)
      , fQ2(fQ2_In)
      , fm2prime(fm2prime_In)
      , fS(fS_In)
      , bD1(bD1_In)
      , bD2(bD2_In)
      , bR(bR_In)
      , bRbar(bRbar_In)
      , bQ1(bQ1_In)
      , bQ2(bQ2_In)
      , bm2prime(bm2prime_In)
      , CollisionMatrix(CollisionMatrix_In)
      , M(std::size(z),
          2 * (nFermions + nBosons),
          2 * (nFermions + nBosons))                    // A^-1 * Gamma matrix
      , STilde(std::size(z), 2 * (nFermions + nBosons)) // A^-1 * Sources vector

  {
    gamma = 1 / sqrt(1 - vw * vw); // Lorenz factor

    for (int k = 0; k < std::size(z); k++)
    {
      // Calculate A^-1
      MatDoub Ainverse(2 * (nBosons + nFermions), 2 * (nBosons + nFermions));
      for (int fermion = 0; fermion < nFermions; fermion++)
      {
        Ainverse[2 * fermion][2 * fermion] =
            fR[k][fermion] /
            (fD2[k][fermion] - fD1[k][fermion] * fR[k][fermion]);
        Ainverse[2 * fermion][2 * fermion + 1] =
            -1 / (fD2[k][fermion] - fD1[k][fermion] * fR[k][fermion]);
        Ainverse[2 * fermion + 1][2 * fermion] =
            fD2[k][fermion] /
            (fD2[k][fermion] - fD1[k][fermion] * fR[k][fermion]);
        Ainverse[2 * fermion + 1][2 * fermion + 1] =
            -fD1[k][fermion] /
            (fD2[k][fermion] - fD1[k][fermion] * fR[k][fermion]);
      }

      for (int boson = 0; boson < nBosons; boson++)
      {
        Ainverse[2 * nFermions + 2 * boson][2 * nFermions + 2 * boson] =
            bR[k][boson] / (bD2[k][boson] - bD1[k][boson] * bR[k][boson]);
        Ainverse[2 * nFermions + 2 * boson][2 * nFermions + 2 * boson + 1] =
            -1 / (bD2[k][boson] - bD1[k][boson] * bR[k][boson]);
        Ainverse[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson] =
            bD2[k][boson] / (bD2[k][boson] - bD1[k][boson] * bR[k][boson]);
        Ainverse[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson + 1] =
            -bD1[k][boson] / (bD2[k][boson] - bD1[k][boson] * bR[k][boson]);
      }

      // Calculate Gamma Matrix (Collision matrix - m^2' B)
      MatDoub Gamma(2 * (nBosons + nFermions), 2 * (nBosons + nFermions));

      // Gamma = CollisionMatrix[k]
      for (int i = 0; i < 2 * (nBosons + nFermions); i++)
        for (int j = 0; j < 2 * (nBosons + nFermions); j++)
          Gamma[i][j] = CollisionMatrix[k][i][j];

      for (int fermion = 0; fermion < nFermions; fermion++)
      {
        Gamma[2 * fermion][2 * fermion] -=
            gamma * vw * fQ1[k][fermion] * fm2prime[k][fermion];
        Gamma[2 * fermion + 1][2 * fermion] -=
            gamma * vw * fQ2[k][fermion] * fm2prime[k][fermion];
        Gamma[2 * fermion + 1][2 * fermion + 1] -=
            fRbar[k][fermion] * fm2prime[k][fermion];
      }
      for (int boson = 0; boson < nBosons; boson++)
      {
        Gamma[2 * nFermions + 2 * boson][2 * nFermions + 2 * boson] -=
            gamma * vw * bQ1[k][boson] * bm2prime[k][boson];
        Gamma[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson] -=
            gamma * vw * bQ2[k][boson] * bm2prime[k][boson];
        ;
        Gamma[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson + 1] -=
            bRbar[k][boson] * bm2prime[k][boson];
      }

      // Calculate M = A^-1 * Gamma
      for (int i = 0; i < 2 * (nBosons + nFermions); i++)
        for (int j = 0; j < 2 * (nBosons + nFermions); j++)
          for (int l = 0; l < 2 * (nBosons + nFermions); l++)
            M[k][i][j] = Ainverse[i][l] * Gamma[l][j];

      // Calculate Stilde = A^-1 * S
      for (int i = 0; i < 2 * (nBosons + nFermions); i++)
        for (int j = 0; j < 2 * (nFermions); j++)
          STilde[k][i] += Ainverse[i][j] * fS[k][j];
    }
  }

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
      for (int j = 0; j < 2 * (nFermions * nBosons); j++)
        for (int n = 0; n < 2 * (nFermions * nBosons); n++)
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
