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

void Q8oint1::set_u(const double u_in)
{
  u = u_in;
}

double Q8oint1::operator()(const double y)
{
  const double w = x + (1 - u) / u;
  return 0.5 * Ki->Kintegrand2D(w, y, x) * f0w(w, s, 1) / (u * u);
}

double Q8oint2::operator()(const double u)
{
  integrand.set_u(u);
  double f[4];
  for (size_t i = 0; i < 4; i++)
    f[i] = integrand(-1. + 2. * (double)i / 3.);
  return adap_simpson38(integrand, -1, 1, f, 1e-4);
}

void Q9oint1::set_u(const double u_in)
{
  u = u_in;
}

double Q9oint1::operator()(const double y)
{
  const double w = x + (1 - u) / u;
  if (part == 1)
    return 0.25 * Ki->Kintegrand2D(w, y, x) * f0w(w, s, 1) / (u * u);
  else
    return 0.25 * Ki->gamw * Ki->Kintegrand2D(w, y, x) * f0w(w, s, 2) / (u * u);
}

double Q9oint2::operator()(const double u)
{
  integrand.set_u(u);
  double f[4];
  for (size_t i = 0; i < 4; i++)
    f[i] = integrand(-1. + 2. * (double)i / 3.);
  return adap_simpson38(integrand, -1, 1, f, 1e-4);
}

double Kfactor::operator()(const K_type type, const double m)
{
  const double x = m / Ki->Tc;
  switch (type)
  {
  case K_type::N0bos:
  {
    N0int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::N0fer:
  {
    N0int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Rbarbos:
  {
    N0int integrand1(Ki, -1, x);
    const double est1 = kronrod_61(integrand1, 0., 1.);
    const double res1 = h_adap_gauss_kronrod_15(integrand1, 0., 1., est1, 1e-4);
    Rbarint integrand2(Ki, -1, x);
    const double est2 = kronrod_61(integrand2, 0., 1.);
    const double res2 = h_adap_gauss_kronrod_15(integrand2, 0., 1., est2, 1e-4);
    return res2 / res1;
  }
  case K_type::Rbarfer:
  {
    N0int integrand1(Ki, 1, x);
    const double est1 = kronrod_61(integrand1, 0., 1.);
    const double res1 = h_adap_gauss_kronrod_15(integrand1, 0., 1., est1, 1e-4);
    Rbarint integrand2(Ki, 1, x);
    const double est2 = kronrod_61(integrand2, 0., 1.);
    const double res2 = h_adap_gauss_kronrod_15(integrand2, 0., 1., est2, 1e-4);
    return res2 / res1;
  }
  case K_type::D0bos:
  {
    D0int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D0fer:
  {
    D0int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D1bos:
  {
    D0int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return -Ki->vw * h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D1fer:
  {
    D0int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return -Ki->vw * h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D2bos:
  {
    D2int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::D2fer:
  {
    D2int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q1bos:
  {
    Q1int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q1fer:
  {
    Q1int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q2bos:
  {
    Q2int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q2fer:
  {
    Q2int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe1bos:
  {
    Qe1int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe1fer:
  {
    Qe1int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe2bos:
  {
    Qe2int integrand(Ki, -1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Qe2fer:
  {
    Qe2int integrand(Ki, 1, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q8o1bos:
  {
    Q8oint2 integrand(Ki, -1, x);
    Ki->set_nmk(-1, 1, 1);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q8o1fer:
  {
    Q8oint2 integrand(Ki, 1, x);
    Ki->set_nmk(-1, 1, 1);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q8o2bos:
  {
    Q8oint2 integrand(Ki, -1, x);
    Ki->set_nmk(0, 2, 1);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q8o2fer:
  {
    Q8oint2 integrand(Ki, 1, x);
    Ki->set_nmk(0, 2, 1);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  case K_type::Q9o1bos:
  {
    double res = 0.;
    Q9oint2 integrand1(Ki, -1, x, 1);
    Q9oint2 integrand2(Ki, -1, x, 2);
    Ki->set_nmk(-1, 3, 1);
    double est = kronrod_61(integrand1, 0., 1.);
    res += h_adap_gauss_kronrod_15(integrand1, 0., 1., est, 1e-4);
    Ki->set_nmk(-1, 2, 2);
    est = kronrod_61(integrand2, 0., 1.);
    res -= h_adap_gauss_kronrod_15(integrand2, 0., 1., est, 1e-4);
    return res;
  }
  case K_type::Q9o1fer:
  {
    double res = 0.;
    Q9oint2 integrand1(Ki, 1, x, 1);
    Q9oint2 integrand2(Ki, 1, x, 2);
    Ki->set_nmk(-1, 3, 1);
    double est = kronrod_61(integrand1, 0., 1.);
    res += h_adap_gauss_kronrod_15(integrand1, 0., 1., est, 1e-4);
    Ki->set_nmk(-1, 2, 2);
    est = kronrod_61(integrand2, 0., 1.);
    res -= h_adap_gauss_kronrod_15(integrand2, 0., 1., est, 1e-4);
    return res;
  }
  case K_type::Q9o2bos:
  {
    double res = 0.;
    Q9oint2 integrand1(Ki, -1, x, 1);
    Q9oint2 integrand2(Ki, -1, x, 2);
    Ki->set_nmk(0, 4, 1);
    double est = kronrod_61(integrand1, 0., 1.);
    res += h_adap_gauss_kronrod_15(integrand1, 0., 1., est, 1e-4);
    Ki->set_nmk(0, 3, 2);
    est = kronrod_61(integrand2, 0., 1.);
    res -= h_adap_gauss_kronrod_15(integrand2, 0., 1., est, 1e-4);
    return res;
  }
  case K_type::Q9o2fer:
  {
    double res = 0.;
    Q9oint2 integrand1(Ki, 1, x, 1);
    Q9oint2 integrand2(Ki, 1, x, 2);
    Ki->set_nmk(0, 4, 1);
    double est1 = kronrod_61(integrand1, 0., 1.);
    Ki->set_nmk(0, 3, 2);
    double est2 = kronrod_61(integrand2, 0., 1.);
    if (Ki->fast) return est1 - est2;
    Ki->set_nmk(0, 4, 1);
    res += h_adap_gauss_kronrod_15(integrand1, 0., 1., est1, 1e-4);
    Ki->set_nmk(0, 3, 2);
    res -= h_adap_gauss_kronrod_15(integrand2, 0., 1., est2, 1e-4);
    return res;
  }

  break;
  default: break;
  }
}

} // namespace BSMPT
