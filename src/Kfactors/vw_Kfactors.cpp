#include <BSMPT/Kfactors/vw_Kfactors.h>

namespace BSMPT
{

double f0w(const double w, const double s, const int diff)
{
  if (w > 100) return std::pow(-1, diff) * std::exp(-w);
  const double expw = std::exp(-w / 2.);
  if (diff == 0)
    return expw / (1 / expw + s * expw);
  else if (diff == 1)
    return -1 / std::pow(1 / expw + s * expw, 2);
  else if (diff == 2)
    return (1 / expw - s * expw) / std::pow(1 / expw + s * expw, 3);
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
  const double Vx  = ((y == -1 || y == 1) && w > 100 * x)
                         ? (1 + Ki->vw) * x * x / ((1 - Ki->vw) * w * w)
                         : 1. / (1. + x * x / (pzt * pzt)) /
                              std::sqrt(1. - pzt * pzt / (Et * Et));
  return Vx / pzt;
}

double Q8o1int2::operator()(const double u)
{
  integrand.set_w(u);
  if (std::abs(integrand.pre) == 0.) return 0.;
  return integrand.pre * adap_gauss_kronrod_15(integrand, -1, 1, 1e-4) /
         (u * u);
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
  const double Vx  = ((y == -1 || y == 1) && w > 100 * x)
                         ? (1 + Ki->vw) * x * x / ((1 - Ki->vw) * w * w)
                         : 1. / (1. + x * x / (pzt * pzt)) /
                              std::sqrt(1. - pzt * pzt / (Et * Et));
  return Vx / Et;
}

double Q8o2int2::operator()(const double u)
{
  integrand.set_w(u);
  if (std::abs(integrand.pre) == 0.) return 0.;
  return integrand.pre * adap_gauss_kronrod_15(integrand, -1, 1, 1e-4) /
         (u * u);
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
  if (std::abs(integrand.pre) == 0.) return 0.;
  return integrand.pre * adap_gauss_kronrod_15(integrand, -1, 1, 1e-8) /
         (u * u);
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
  if (std::abs(integrand.pre) == 0.) return 0.;
  return integrand.pre * adap_gauss_kronrod_15(integrand, -1, 1, 1e-4) /
         (u * u);
}

void Kfactor::spline_Kfactors(const double mmin,
                              const double mmax,
                              const size_t N_points)
{
  double logmmin = std::log10(mmin);
  double logmmax = std::log10(mmax);
  double h       = (logmmax - logmmin) / (double)N_points;
  std::vector<double> masses;
  std::vector<double> K_values;
  for (size_t i = 0; i < 14; i++)
  {
    masses.clear();
    K_values.clear();
    for (double mass = logmmin; mass <= logmmax; mass += h)
    {
      double current = pow(10, mass);
      masses.push_back(current);
      K_values.push_back(operator()((K_type)i, fermion, current));
    }
    Kspl[i] = tk::spline(masses, K_values);
  }

  // only until Qe1 in the boson case since there are no source terms
  for (size_t i = 0; i < 9; i++)
  {
    std::vector<double> masses;
    std::vector<double> K_values;
    for (double mass = logmmin; mass <= logmmax; mass += h)
    {
      double current = pow(10, mass);
      masses.push_back(current);
      K_values.push_back(operator()((K_type)i, boson, current));
    }
    Kspl[14 + i] = tk::spline(masses, K_values);
  }
  is_spline = true;
}

double
Kfactor::operator()(const K_type ktype, const P_type ptype, const double m)
{
  const double x      = m / Ki->Tc;
  const int statistic = (ptype == boson ? -1 : 1);
  const size_t shift  = (ptype == boson ? 14 : 0);
  switch (ktype)
  {
  case K_type::N0:
  {
    if (is_spline) return Kspl[shift + 0](m);
    N0int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Rbar:
  {
    if (is_spline) return Kspl[shift + 1](m);
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
    if (is_spline) return Kspl[shift + 2](m);
    K4int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D0:
  {
    if (is_spline) return Kspl[shift + 3](m);
    D0int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D1:
  {
    if (is_spline) return Kspl[shift + 4](m);
    D0int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return -Ki->vw * h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D2:
  {
    if (is_spline) return Kspl[shift + 5](m);
    D2int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q1:
  {
    if (is_spline) return Kspl[shift + 6](m);
    Q1int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q2:
  {
    if (is_spline) return Kspl[shift + 7](m);
    Q2int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe1:
  {
    if (is_spline) return Kspl[shift + 8](m);
    Qe1int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe2:
  {
    if (is_spline) return Kspl[shift + 9](m);
    Qe2int integrand(Ki, statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q8o1:
  {
    if (is_spline) return Kspl[shift + 10](m);
    Q8o1int2 integrand(Ki, statistic, x);
    const double res = adap_gauss_kronrod_15(integrand, 0., 1., 1e-4);
    return -3 / (2. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -2) * res;
  }
  case K_type::Q8o2:
  {
    if (is_spline) return Kspl[shift + 11](m);
    Q8o2int2 integrand(Ki, statistic, x);
    const double res = adap_gauss_kronrod_15(integrand, 0., 1., 1e-4);
    return -3 / (2. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -2) * res;
  }
  case K_type::Q9o1:
  {
    if (is_spline) return Kspl[shift + 12](m);
    Q9o1int2 integrand1(Ki, statistic, x, 1);
    const double res1 = adap_gauss_kronrod_15(integrand1, 0., 1., 1e-4);
    Q9o1int2 integrand2(Ki, statistic, x, 2);
    const double res2 = adap_gauss_kronrod_15(integrand2, 0., 1., 1e-4);

    return -3 / (4. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -4) * (res1 - res2);
  }
  case K_type::Q9o2:
  {
    if (is_spline) return Kspl[shift + 13](m);
    Q9o2int2 integrand1(Ki, statistic, x, 1);
    Q9o2int2 integrand2(Ki, statistic, x, 2);
    const double res1 = adap_gauss_kronrod_15(integrand1, 0., 1., 1e-4);
    const double res2 = adap_gauss_kronrod_15(integrand2, 0., 1., 1e-4);

    return -3 / (4. * M_PI * M_PI * Ki->gamw) * pow(Ki->Tc, -4) * (res1 - res2);
  }

  break;
  default: break;
  }
}

} // namespace BSMPT
