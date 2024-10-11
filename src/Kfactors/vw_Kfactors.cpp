#include <BSMPT/Kfactors/vw_Kfactors.h>

namespace BSMPT
{

double f0w(const double w, const double s, const int diff)
{
  if (w > 100) return std::pow(-1, diff) * std::exp(-w);
  double res        = 0.;
  const double expw = std::exp(-w / 2.);
  if (diff == 0)
    res += expw / (1 / expw + s * expw);
  else if (diff == 1)
    res -= 1 / std::pow(1 / expw + s * expw, 2);
  else if (diff == 2)
    res += (1 / expw - s * expw) / std::pow(1 / expw + s * expw, 3);
  return res;
}

double N0int::operator()(const double u)
{
  const double w = x + (1 - u) / u;
  return 4 * M_PI * Ki->Tc * Ki->Tc * Ki->Tc * std::sqrt(w * w - x * x) * w *
         f0w(w, s, 0) / (u * u);
}

double Rbarint::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double pwt = std::sqrt(w * w - x * x);
  return Ki->Tc * M_PI / (Ki->gamw * Ki->gamw) *
         log(std::abs((pwt - Ki->vw * w) / (pwt + Ki->vw * w))) * f0w(w, s, 0) /
         (u * u);
}

double K4int::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double pwt = std::sqrt(w * w - x * x);
  return -2 / (M_PI * M_PI) * pwt * pwt * pwt / w * f0w(w, s, 1) / (u * u);
}

double D0int::operator()(const double u)
{
  const double w = x + (1 - u) / u;
  return -6. / (M_PI * M_PI) * w * std::sqrt(w * w - x * x) * f0w(w, s, 1) /
         (u * u);
}

double D2int::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double vw2 = Ki->vw * Ki->vw;
  const double pwt = std::sqrt(w * w - x * x);
  return -6 / (M_PI * M_PI * Ki->vw * vw2) * w *
         (Ki->vw * pwt * (2 * vw2 - 1) +
          (vw2 - 1) * (vw2 - 1) * w * std::atanh(Ki->vw * pwt / w)) *
         f0w(w, s, 1) / (u * u);
}

double Q1int::operator()(const double u)
{
  const double w = x + (1 - u) / u;
  return -3 / (Ki->gamw * M_PI * M_PI * Ki->Tc * Ki->Tc) *
         std::sqrt(w * w - x * x) * f0w(w, s, 2) / (u * u);
}

double Q2int::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double vw2 = Ki->vw * Ki->vw;
  const double pwt = std::sqrt(w * w - x * x);
  return 3 / (Ki->gamw * M_PI * M_PI * Ki->Tc * Ki->Tc * vw2) *
         (Ki->vw * pwt + (vw2 - 1) * w * atanh(Ki->vw * pwt / w)) *
         f0w(w, s, 2) / (u * u);
}

double Qe1int::operator()(const double u)
{
  const double w = x + (1 - u) / u;
  return -3 / (Ki->gamw * M_PI * M_PI * Ki->Tc) * std::sqrt(w * w - x * x) *
         f0w(w, s, 1) / (u * u);
}

double Qe2int::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double vw2 = Ki->vw * Ki->vw;
  const double pwt = std::sqrt(w * w - x * x);
  return 3. / (Ki->gamw * M_PI * M_PI * Ki->Tc * vw2) *
         (Ki->vw * pwt + (vw2 - 1) * w * atanh(Ki->vw * pwt / w)) *
         f0w(w, s, 1) / (u * u);
}

void Q8o1int1::set_w(const double u_in)
{
  w   = x + (1 - u_in) / u_in;
  pwt = std::sqrt(w * w - x * x);
  pre = pwt * f0w(w, s, 1);
}

double Q8o1int1::operator()(const double y)
{
  const double pzt = Ki->gamw * (y * pwt - w * Ki->vw);
  const double Et  = Ki->gamw * (w - Ki->vw * y * pwt);
  const double Vx =
      1. / (1. + x * x / (pzt * pzt)) / sqrt(1. - pzt * pzt / (Et * Et));
  return Vx / pzt;
}

double Q8o1int2::operator()(const double u)
{
  integrand.set_w(u);
  double f[4];
  for (size_t i = 0; i < 4; i++)
    f[i] = integrand(-1. + 2. * (double)i / 3.);
  return integrand.pre * adap_simpson38(integrand, -1, 1, f, 1e-4) / (u * u);
}

void Q8o2int1::set_w(const double u_in)
{
  w   = x + (1 - u_in) / u_in;
  pwt = std::sqrt(w * w - x * x);
  pre = pwt * f0w(w, s, 1);
}

double Q8o2int1::operator()(const double y)
{
  const double pzt = Ki->gamw * (y * pwt - w * Ki->vw);
  const double Et  = Ki->gamw * (w - Ki->vw * y * pwt);
  const double Vx =
      1. / (1. + x * x / (pzt * pzt)) / sqrt(1. - pzt * pzt / (Et * Et));
  return Vx / Et;
}

double Q8o2int2::operator()(const double u)
{
  integrand.set_w(u);
  double f[4];
  for (size_t i = 0; i < 4; i++)
    f[i] = integrand(-1. + 2. * (double)i / 3.);
  return adap_simpson38(integrand, -1, 1, f, 1e-4) * integrand.pre / (u * u);
}

void Q9o1int1::set_w(const double u_in)
{
  w   = x + (1 - u_in) / u_in;
  pwt = std::sqrt(w * w - x * x);
  if (part == 1)
    pre = f0w(w, s, 1) * pwt;
  else
    pre = Ki->gamw * f0w(w, s, 2) * pwt;
}

double Q9o1int1::operator()(const double y)
{
  const double pzt = Ki->gamw * (y * pwt - w * Ki->vw);
  const double Et  = Ki->gamw * (w - Ki->vw * y * pwt);
  const double Vx =
      1. / (1. + x * x / (pzt * pzt)) / sqrt(1. - pzt * pzt / (Et * Et));
  if (part == 1)
    return Vx / (pzt * Et * Et);
  else
    return Vx / (pzt * Et);
}

double Q9o1int2::operator()(const double u)
{
  integrand.set_w(u);
  double f[4];
  for (size_t i = 0; i < 4; i++)
    f[i] = integrand(-1. + 2. * (double)i / 3.);
  return integrand.pre * adap_simpson38(integrand, -1, 1, f, 1e-4) / (u * u);
}

void Q9o2int1::set_w(const double u_in)
{
  w   = x + (1 - u_in) / u_in;
  pwt = std::sqrt(w * w - x * x);
  if (part == 1)
    pre = pwt * f0w(w, s, 1);
  else
    pre = pwt * Ki->gamw * f0w(w, s, 2);
}

double Q9o2int1::operator()(const double y)
{
  const double pzt = Ki->gamw * (y * pwt - w * Ki->vw);
  const double Et  = Ki->gamw * (w - Ki->vw * y * pwt);
  const double Vx =
      1. / (1. + x * x / (pzt * pzt)) / sqrt(1. - pzt * pzt / (Et * Et));
  if (part == 1)
    return Vx / (Et * Et * Et);
  else
    return Vx / (Et * Et);
}

double Q9o2int2::operator()(const double u)
{
  integrand.set_w(u);
  double f[4];
  for (size_t i = 0; i < 4; i++)
    f[i] = integrand(-1. + 2. * (double)i / 3.);
  return integrand.pre * adap_simpson38(integrand, -1, 1, f, 1e-4) / (u * u);
}

double
Kfactor::operator()(const K_type ktype, const P_type ptype, const double m)
{
  const double x      = m / Ki->Tc;
  const int statistic = (ptype == boson ? -1 : 1);
  switch (ktype)
  {
  case K_type::N0:
  {
    N0int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Rbar:
  {
    N0int integrand1(Ki, statistic, x);
    const double est1 = kronrod_61(integrand1, 0., 1.);
    const double res1 = h_adap_gauss_kronrod_15(integrand1, 0., 1., est1, 1e-4);
    Rbarint integrand2(Ki, statistic, x);
    const double est2 = kronrod_61(integrand2, 0., 1.);
    const double res2 = h_adap_gauss_kronrod_15(integrand2, 0., 1., est2, 1e-4);
    return res2 / res1;
  }
  case K_type::K4:
  {
    K4int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D0:
  {
    D0int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D1:
  {
    D0int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return -Ki->vw * h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D2:
  {
    D2int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q1:
  {
    Q1int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q2:
  {
    Q2int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe1:
  {
    Qe1int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe2:
  {
    Qe2int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q8o1:
  {
    Q8o1int2 integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    const double res = h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
    return -3 / (2. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -2) * res;
  }
  case K_type::Q8o2:
  {
    Q8o2int2 integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    const double res = h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
    return -3 / (2. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -2) * res;
  }
  case K_type::Q9o1:
  {
    Q9o1int2 integrand1(Ki, statistic, x, 1);
    const double est1 = kronrod_61(integrand1, 0., 1.);
    const double res1 = h_adap_gauss_kronrod_15(integrand1, 0., 1., est1, 1e-4);
    Q9o1int2 integrand2(Ki, statistic, x, 2);
    const double est2 = kronrod_61(integrand2, 0., 1.);
    const double res2 = h_adap_gauss_kronrod_15(integrand2, 0., 1., est2, 1e-4);

    return -3 / (4. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -4) * (res1 - res2);
  }
  case K_type::Q9o2:
  {
    Q9o2int2 integrand1(Ki, statistic, x, 1);
    const double est1 = kronrod_61(integrand1, 0., 1.);
    const double res1 = h_adap_gauss_kronrod_15(integrand1, 0., 1., est1, 1e-4);
    Q9o2int2 integrand2(Ki, statistic, x, 2);
    const double est2 = kronrod_61(integrand2, 0., 1.);
    const double res2 = h_adap_gauss_kronrod_15(integrand2, 0., 1., est2, 1e-4);

    return -3 / (4. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -4) * (res1 - res2);
  }

  break;
  default: break;
  }
}

} // namespace BSMPT
