#include <BSMPT/baryo_calculation/transport_network.h>

namespace BSMPT
{

P_type TransportNetwork::get_particle_type(const Particles prtcl)
{
  if ((prtcl == tL) || (prtcl == bL)) return fermion;
  if (prtcl == tR) return antifermion;
  return boson;
}

double TransportNetwork::vev1_profile(const double &z, const size_t deriv)
{
  double th = tanh(z / LW);
  if (deriv == 0)
    return 0.5 * (1 - th);
  else if (deriv == 1)
    return 0.5 / LW * (th * th - 1);
  else if (deriv == 2)
    return th / (LW * LW) * (1 - th * th);
  else
    return 0.;
}

double TransportNetwork::top_mass(const double z, const size_t deriv)
{
  double th    = tanh(z / LW);
  double hz0   = vh / 2. * (1. - th);
  double sz0   = vs / 2. * (1. + th);
  double phase = std::sqrt(1. + sz0 * sz0 / (LAM * LAM));
  if (deriv == 0)
    return 0.7 * hz0 * phase;
  else if (deriv == 1)
  {
    double hz1 = vh / (2. * LW) * (th * th - 1);
    double sz1 = vs / (2. * LW) * (1 - th * th);
    return 0.7 * (hz1 * phase + hz0 * sz0 * sz1 / (phase * LAM * LAM));
  }
  else
    return 0.;
}

double TransportNetwork::W_mass(const double z)
{
  static const double g = 0.641;
  double hz0            = vh / 2. * (1. - tanh(z / LW));
  return g * hz0 / 2.;
}

double TransportNetwork::theta(const double z, const size_t deriv)
{
  double th  = tanh(z / LW);
  double sz0 = vs / 2. * (1. + th);
  if (deriv == 0)
    return atan(sz0 / LAM);
  else if (deriv == 1)
  {
    double sz1 = vs / (2. * LW) * (1 - th * th);
    return LAM * sz1 / (LAM * LAM + sz0 * sz0);
  }
  else if (deriv == 2)
  {
    double sz1 = vs / (2. * LW) * (1 - th * th);
    double sz2 = vs / (LW * LW) * th * (th * th - 1.);
    return LAM * ((LAM * LAM + sz0 * sz0) * sz2 - 2 * sz0 * sz1 * sz1) /
           ((LAM * LAM + sz0 * sz0) * (LAM * LAM + sz0 * sz0));
  }
}

void TransportNetwork::spline_Kfactors(const double zmin,
                                       const double zmax,
                                       const size_t N_points)
{
  double mmax = top_mass(zmin, 0);
  double mmin = top_mass(zmax, 0);
  K.spline_Kfactors(mmin, mmax, N_points);
}

VecDoub TransportNetwork::calc_Source(const double z)
{
  VecDoub res(2);
  const double dth  = theta(z, 1);
  const double d2th = theta(z, 2);
  double mt         = top_mass(z, 0);
  double dmt        = top_mass(z, 1);
  double p1         = (2 * mt * dmt * dth + mt * mt * d2th);
  double p2         = 2 * mt * dmt * mt * mt * dth;
  res[0]            = Ki->vw * Ki->gamw *
           (p1 * K(Q8o1, fermion, mt) - p2 * K(Q9o1, fermion, mt));
  res[1] = Ki->vw * Ki->gamw *
           (p1 * K(Q8o2, fermion, mt) - p2 * K(Q9o2, fermion, mt));
  return res;
}

void TransportNetwork::operator()(const double z, VecDoub &u, VecDoub &du)
{
  double mt  = top_mass(z, 0);
  double dmt = top_mass(z, 1);

  double D0t = K(D0, fermion, mt);
  double D1t = K(D1, fermion, mt);
  double D2t = K(D2, fermion, mt);
  double Q1t = K(Q1, fermion, mt);
  double Q2t = K(Q2, fermion, mt);
  double Rbt = K(Rbar, fermion, mt);

  double detAt = Ki->vw * D1t + D2t;

  double S1, S2;
  {
    VecDoub temp = calc_Source(z);
    S1           = temp[0];
    S2           = temp[1];
  }
  double Gh    = pow(W_mass(z), 2) / (50. * Ki->Tc) * u[6];
  double Gtott = K(K4, fermion, mt) * Ki->Tc / (6. * D0t);
  double GW    = Gtoth * (u[2] - u[0]);
  double GY    = 4.2e-3 * Ki->Tc;
  double GM    = mt * mt / (63. * Ki->Tc) * (u[4] - u[0]);
  double GSS   = 4.9e-4 * Ki->Tc *
               ((1. + 9. * D0t) * u[0] + (1. + 9. * D0b) * u[2] -
                (1. - 9. * D0t) * u[4]);

  du[0] = (2. * mt * dmt *
               (Rbt * u[1] + Ki->gamw * Ki->vw * (Ki->vw * Q1t + Q2t) * u[0]) +
           Ki->vw * (GW + GM + GY * (u[4] - u[6] - u[0]) - GSS) + Gtott * u[1] -
           Ki->vw * S1 - S2) /
          detAt;
  du[1] = (2. * mt * dmt *
               (Rbt * D1t * u[1] +
                Ki->gamw * Ki->vw * (D1t * Q2t - D2t * Q1t) * u[0]) -
           D2t * (GW + GM + GY * (u[4] - u[6] - u[0]) - GSS) +
           D1t * Gtott * u[1] + S1 * D2t - S2 * D1t) /
          detAt;
  du[2] =
      (Ki->vw * (-GW + GY * (u[4] - u[6] - u[2]) - GSS) + Gtotb * u[3]) / detAb;
  du[3] = (D2b * (GW + GY * (u[2] + u[6] - u[4]) + GSS) + D1b * Gtotb * u[3]) /
          detAb;
  du[4] = (2 * mt * dmt *
               (Rbt * u[5] + Ki->gamw * Ki->vw * (Ki->vw * Q1t + Q2t) * u[4]) +
           Ki->vw * (-GM + GY * (u[2] + 2. * u[6] + u[0] - 2. * u[4]) + GSS) +
           Gtoth * u[5] + Ki->vw * S1 + S2) /
          detAt;
  du[5] = (2. * mt * dmt *
               (Rbt * D1t * u[5] +
                Ki->gamw * Ki->vw * (D1t * Q2t - D2t * Q1t) * u[4]) +
           D2t * (GM - GY * (u[2] + 2. * u[6] + u[0] - 2. * u[4]) - GSS) +
           D1t * Gtott * u[5] + S2 * D1t - S1 * D2t) /
          detAt;
  du[6] = (-Ki->vw * (Gh + GY * (u[2] + 2. * u[6] + u[0] - 2. * u[4])) +
           Gtoth * u[7]) /
          detAh;
  du[7] = (D2h * (Gh + GY * (u[2] + 2. * u[6] + u[0] - 2. * u[4])) +
           D1h * Gtoth * u[7]) /
          detAh;
}

VecDoub shootf::operator()(VecDoub &v)
{
  VecDoub ul = {0., v[0], 0., v[1], 0., v[2], 0., v[3]};
  VecDoub ur = {0., v[4], 0., v[5], 0., v[6], 0., v[7]};
  VecDoub dul(8), dur(8);
  VecDoub res(8);
  double zini = zl;
  double zfin = zr;
  rk4_adap(tr, zini, ul, zmid, 1e-4, 1e-4, 1e-5, save);
  tr(zmid, ul, dul);
  // if (save) rk4_adap(tr, zini, ul, zfin, 1e-4, 1e-4, 1e-4, save);
  rk4_adap(tr, zfin, ur, zmid, -1e-4, 1e-4, -1e-5, save);
  tr(zmid, ur, dur);
  for (size_t i = 0; i < 4; i++)
  {
    res[2 * i]     = 2*(ul[2 * i] - ur[2 * i]);
    res[2 * i + 1] = dul[2 * i] + dur[2 * i];
    if (save)
      std::cout << "i: " << i << " " << dul[2 * i] << " " << dur[2 * i] << "\n";
  }
  if (save)
    for (auto it : res)
      std::cout << it << "\n";
  // if (save) rk4_adap(tr, zfin, ur, zini, -1e-4, 1e-4, -1e-4, save);

  return res;
}

/* VecDoub shootf::operator()(VecDoub &v)
{
  VecDoub ur = {v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]};
  VecDoub dur(8);
  VecDoub res(8);
  double zini = zl;
  double zfin = zr;
  tr(zfin, ur, dur);
  rk4_adap(tr, zfin, ur, 0., -1e-4, 1e-4, -1e-5, save);
  for (size_t i = 0; i < 4; i++)
  {
    res[2 * i]     = std::abs(ur[2 * i]) + std::abs(dur[2 * i]);
    res[2 * i + 1] = std::abs(dur[2 * i + 1]);
  }
  if (save)
    for (auto it : res)
      std::cout << it << "\n";
  // if (save) rk4_adap(tr, zfin, ur, zini, -1e-4, 1e-4, -1e-4, save);

  return res;
} */
} // namespace BSMPT