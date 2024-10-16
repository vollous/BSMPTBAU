#include <BSMPT/baryo_calculation/difeq.h>

namespace BSMPT
{
Difeq::Difeq(
    // z list
    VecDoub &z_In,
    TransportNetwork tn_in)
    : z(z_In)
    , TN(tn_in)
{
  double zmid;
  MatDoub Ainv;
  MatDoub Gamma;
  VecDoub Source;
  for (size_t i = 1; i < z.size(); i++)
  {
    zmid = 0.5 * (z[i] + z[i - 1]);

    // Ainv * (Collision - B)
    Ainv  = TN.calc_A_inv(zmid);
    Gamma = TN.calc_Collision(zmid) - TN.calc_B(zmid);
    M.push_back(Ainv * Gamma);

    // Ainv * Source
    Source = TN.calc_Source(zmid);
    Stilde.push_back(Ainv * Source);
  }
}

void Difeq::smatrix(const int k,
                    const int k1,
                    const int k2,
                    const int jsf,
                    const int is1,
                    const int isf,
                    VecInt &indexv,
                    MatDoub &s,
                    MatDoub &y)
{
  size_t Np = TN.get_N_particles();
  if (k == k1)
  {
    set_zero(s);
    // Boundary conditions mu = 0 on first boundary
    for (size_t i = 0; i < Np; i++)
    {
      // Sn at the first boundary
      s[Np + i][2 * Np + indexv[2 * i]] = 1.0;
      // B0
      s[Np + i][jsf] = y[indexv[2 * i]][0];
    }
  }
  else if (k > k2 - 1)
  {
    set_zero(s);
    // Boundary conditions mu = 0 on second boundary
    for (size_t i = 0; i < Np; i++)
    {
      // Sn at the last boundary
      s[i][2 * Np + indexv[2 * i]] = 1.0;
      // C0
      s[i][jsf] = y[indexv[2 * i]][z.size()];
    }
  }
  else
  {
    double temp;
    for (size_t j = 0; j < 2 * Np; j++)
      for (size_t n = 0; n < 2 * Np; n++)
      {
        double del = (j == n) ? 1 : 0;
        // s matrix for the middle point
        // S_{j,n}
        s[j][indexv[n]] = -del - 0.5 * (z[k] - z[k - 1]) * (M[k - 1][j][n]);
        // S_{j,N + n}
        s[j][2 * Np + indexv[n]] =
            del - 0.5 * (z[k] - z[k - 1]) * (M[k - 1][j][n]);
        // Equations for E(k,k-1)
        temp = Stilde[k - 1][j];
        for (size_t i = 0; i < 2 * Np; i++)
          temp += M[k - 1][j][i] * (y[i][k] + y[i][k - 1]);
        s[j][jsf] = y[j][k] - y[j][k - 1] - 0.5 * (z[k] - z[k - 1]) * temp;
      }
  }
}

} // namespace BSMPT
