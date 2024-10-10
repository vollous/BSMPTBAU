#include <BSMPT/baryo_calculation/difeq.h>

namespace BSMPT
{
Difeq::Difeq(std::shared_ptr<Kinfo> K_in,
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
             const Mat3DDoub &CollisionMatrix)
    : z(z_In)
    , nFermions(nFermions_In)
    , nBosons(nBosons_In)
    , M(z.size(),
        2 * (nFermions + nBosons),
        2 * (nFermions + nBosons))                // A^-1 * Gamma matrix
    , STilde(z.size(), 2 * (nFermions + nBosons)) // A^-1 * Sources vector

{
  Ki = K_in;
  Kfactor K(Ki);
  double zmid;
  for (int k = 1; k < z.size(); k++)
  {
    double zm = 0.5 * (z[k] + z[k - 1]);
    // Calculate A^-1
    MatDoub Ainverse(2 * (nBosons + nFermions), 2 * (nBosons + nFermions));
    for (int fermion = 0; fermion < nFermions; fermion++)
    {
      Ainverse[2 * fermion][2 * fermion] =
          -Ki->vw / (fD2[k][fermion] + Ki->vw * fD1[k][fermion]);
      Ainverse[2 * fermion][2 * fermion + 1] =
          -1 / (fD2[k][fermion] + Ki->vw * fD1[k][fermion]);
      Ainverse[2 * fermion + 1][2 * fermion] =
          fD2[k][fermion] / (fD2[k][fermion] + Ki->vw * fD1[k][fermion]);
      Ainverse[2 * fermion + 1][2 * fermion + 1] =
          -fD1[k][fermion] / (fD2[k][fermion] + Ki->vw * fD1[k][fermion]);
    }

    for (int boson = 0; boson < nBosons; boson++)
    {
      Ainverse[2 * nFermions + 2 * boson][2 * nFermions + 2 * boson] =
          -Ki->vw / (bD2[k][boson] + Ki->vw * bD1[k][boson]);
      Ainverse[2 * nFermions + 2 * boson][2 * nFermions + 2 * boson + 1] =
          -1 / (bD2[k][boson] + Ki->vw * bD1[k][boson]);
      Ainverse[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson] =
          bD2[k][boson] / (bD2[k][boson] + Ki->vw * bD1[k][boson]);
      Ainverse[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson + 1] =
          -bD1[k][boson] / (bD2[k][boson] + Ki->vw * bD1[k][boson]);
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
          Ki->gamw * Ki->vw * fQ1[k][fermion] * fm2prime[k][fermion];
      Gamma[2 * fermion + 1][2 * fermion] -=
          Ki->gamw * Ki->vw * fQ2[k][fermion] * fm2prime[k][fermion];
      Gamma[2 * fermion + 1][2 * fermion + 1] -=
          fRbar[k][fermion] * fm2prime[k][fermion];
    }
    for (int boson = 0; boson < nBosons; boson++)
    {
      Gamma[2 * nFermions + 2 * boson][2 * nFermions + 2 * boson] -=
          Ki->gamw * Ki->vw * bQ1[k][boson] * bm2prime[k][boson];
      Gamma[2 * nFermions + 2 * boson + 1][2 * nFermions + 2 * boson] -=
          Ki->gamw * Ki->vw * bQ2[k][boson] * bm2prime[k][boson];
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

} // namespace BSMPT
