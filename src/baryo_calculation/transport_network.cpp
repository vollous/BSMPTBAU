#include <BSMPT/baryo_calculation/transport_network.h>

namespace BSMPT
{

size_t TransportNetwork::get_N_particles()
{
  return prtcl_list.size();
}

P_type TransportNetwork::get_particle_type(const Particles prtcl)
{
  if ((prtcl == tL) || (prtcl == bL)) return fermion;
  if (prtcl == tR) return antifermion;
  return boson;
}

double TransportNetwork::vevProfileKink(const double &z, size_t deriv)
{
  double th = tanh(z / LW);

  double res = 0;
  if (deriv == 0)
    res += 0.5 * (1 - th);
  else if (deriv == 1)
    res += 0.5 / LW * (th * th - 1);
  else if (deriv == 2)
    res += th / (LW * LW) * (1 - th * th);
  return res;
}

VecDoub TransportNetwork::calc_vev(const double &z, size_t deriv)
{
  VecDoub vev;
  double val = vevProfileKink(z, deriv);

  for (std::size_t i = 0; i < modelPointer->get_NHiggs(); i++)
  {
    vev.push_back(vev_critical.at(i) * val);
  }
  return vev;
}

VecDoub TransportNetwork::get_top_mass_and_derivative(const VecDoub &vev) const
{
  VecDoub restmp;
  VecDoub res;
  int postop            = modelPointer->get_NQuarks() - 1;
  int posdifftop        = 2 * modelPointer->get_NQuarks() - 1;
  double topmasssquared = -1;
  VecDoub topderivatives;
  for (std::size_t i = 1; i <= modelPointer->get_NHiggs(); i++)
  {
    restmp.clear();
    restmp = modelPointer->QuarkMassesSquared(vev, i);
    if (topmasssquared == -1) topmasssquared = restmp.at(postop);
    topderivatives.push_back(restmp.at(posdifftop));
  }
  res.push_back(topmasssquared);
  if (std::isnan(topmasssquared))
  {
    std::stringstream ss;
    for (std::size_t i = 0; i < vev.size(); i++)
      ss << vev.at(i) << "\t";
    Logger::Write(LoggingLevel::EWBGDetailed, ss.str());
    std::string retmessage = "Nan found in ";
    retmessage += __func__;
    throw std::runtime_error(retmessage);
  }
  for (std::size_t i = 0; i < topderivatives.size(); i++)
  {
    res.push_back(topderivatives.at(i));
    if (std::isnan(topderivatives.at(i)))
    {
      std::string retmessage = "Nan found in ";
      retmessage += __func__;
      retmessage += " at deriv number ";
      retmessage += std::to_string(i);
      std::stringstream ss;
      for (std::size_t j = 0; j < vev.size(); j++)
        ss << vev.at(j) << "\t";
      Logger::Write(LoggingLevel::EWBGDetailed, ss.str());
      throw std::runtime_error(retmessage);
    }
  }
  return res;
}

VecDoub TransportNetwork::get_h_mass_and_derivative(const VecDoub &vev)
{
  VecDoub res;
  res.push_back(modelPointer->HiggsMassesSquared(vev, Ki->Tc).at(0));

  for (size_t i = 1; i <= vev.size(); i++)
  {
    res.push_back(modelPointer->HiggsMassesSquared(vev, Ki->Tc, i).at(0));
  }

  return res;
}

double TransportNetwork::get_W_mass(const VecDoub &vev) const
{
  VecDoub res;
  res = modelPointer->GaugeMassesSquared(vev, Ki->Tc);
  VecDoub nrepeat(modelPointer->get_NGauge());
  for (std::size_t i = 0; i < modelPointer->get_NGauge(); i++)
  {
    nrepeat[i] = 0;
    for (std::size_t j = 0; j < modelPointer->get_NGauge(); j++)
    {
      if (std::abs(res.at(i) - res.at(j)) <= 1e-5) nrepeat[i]++;
    }
  }

  for (int j = modelPointer->get_NGauge() - 1; j >= 0; j--)
  {
    if (nrepeat[j] > 1)
    {
      if (std::isnan(res.at(j)))
      {
        std::string retmessage = "Nan found in ";
        retmessage += __func__;
        throw std::runtime_error(retmessage);
      }
      return res.at(j);
    }
  }
  return 0;
}

double TransportNetwork::calculate_theta(const double &z, size_t diff)
{
  double res = 0;
  // double thetasym = symmetric_CP_violating_phase;
  double thetabrk  = thbrk;
  double difftheta = thbrk - thsym;
  double vevp      = vevProfileKink(z, diff);
  if (diff == 0)
  {
    res = thetabrk - difftheta * vevp;
  }
  else
  {
    res = difftheta * vevp;
  }

  if (std::isnan(res))
  {
    std::string retmessage = "Nan found in ";
    retmessage += __func__;
    throw std::runtime_error(retmessage);
  }

  // double v3c = vev_critical.at(7), v2c = vev_critical.at(6),
  //        v1c = vev_critical.at(4);
  // double tanbeta_temp_sq =
  //     (std::pow(v3c, 2) + std::pow(v2c, 2)) / std::pow(v1c, 2);
  //
  // if (UseTanBetaSuppression)
  //   res *= 1.0 / (1 + tanbeta_temp_sq); // Eq 27 of 0605242

  return res;
}

VecDoub TransportNetwork::get_squared_mass_and_deriv(const double z,
                                                     const Particles prtcl)
{
  VecDoub res = {0., 0.};
  VecDoub vev = calc_vev(z, 0);
  switch (prtcl)
  {
  case Particles::tL:
  case Particles::tR:
  {
    VecDoub vevdiff = calc_vev(z, 1);
    VecDoub topres  = get_top_mass_and_derivative(vev);
    res[0]          = topres[0];
    for (size_t i = 0; i < vevdiff.size(); i++)
    {
      res[1] += vevdiff[i] * topres[i + 1];
    }
  }
  break;
  case Particles::bL:
  case Particles::bR:
  {
    VecDoub vevdiff = calc_vev(z, 1);
    VecDoub topres  = get_top_mass_and_derivative(vev);
    res[0]          = 0.000581943 * topres[0];
    for (size_t i = 0; i < vevdiff.size(); i++)
    {
      res[1] += 0.000581943 * vevdiff[i] * topres[i + 1];
    }
  }
  break;
  case Particles::h:
  {
    VecDoub vevdiff = calc_vev(z, 1);
    VecDoub hres    = get_h_mass_and_derivative(vev);
    res[0]          = hres[0];
    for (size_t i = 0; i < vevdiff.size(); i++)
    {
      res[1] += vevdiff[i] * hres[i + 1];
    }
  }
  break;
  case Particles::W:
  {
    res[0] = get_W_mass(vev);
  }
  break;

  default: break;
  }
  return res;
}

MatDoub TransportNetwork::calc_A_inv(const double z)
{
  size_t N_prtcls = prtcl_list.size();
  MatDoub res(2 * N_prtcls, VecDoub(2 * N_prtcls, 0));
  double m, detA, D1, D2;
  P_type pt;
  Kfactor K(Ki, false);

  for (size_t i = 0; i < N_prtcls; i++)
  {
    m    = std::sqrt(get_squared_mass_and_deriv(z, prtcl_list[i])[0]);
    pt   = get_particle_type(prtcl_list[i]);
    D1   = K(K_type::D1, pt, m);
    D2   = K(K_type::D2, pt, m);
    detA = D1 + Ki->vw * D2;
    res[2 * i][2 * i]         = -Ki->vw / detA;
    res[2 * i][2 * i + 1]     = -1 / detA;
    res[2 * i + 1][2 * i]     = D2 / detA;
    res[2 * i + 1][2 * i + 1] = -D1 / detA;
  }
  return res;
}

MatDoub TransportNetwork::calc_B(const double z)
{
  size_t N_prtcls = prtcl_list.size();
  MatDoub res(2 * N_prtcls, VecDoub(2 * N_prtcls, 0));
  double m, dmsq;
  VecDoub temp;
  P_type pt;
  Kfactor K(Ki, false);

  for (size_t i = 0; i < N_prtcls; i++)
  {
    temp                      = get_squared_mass_and_deriv(z, prtcl_list[i]);
    m                         = std::sqrt(temp[0]);
    dmsq                      = temp[1];
    pt                        = get_particle_type(prtcl_list[i]);
    res[2 * i][2 * i]         = dmsq * Ki->vw * Ki->gamw * K(K_type::Q1, pt, m);
    res[2 * i + 1][2 * i]     = dmsq * Ki->vw * Ki->gamw * K(K_type::Q2, pt, m);
    res[2 * i + 1][2 * i + 1] = dmsq * K(K_type::Rbar, pt, m);
  }
  return res;
}

MatDoub TransportNetwork::calc_Collision(const double z)
{
  MatDoub res;
  Kfactor K(Ki, false);
  double K0         = 1.;
  const double mtsq = get_squared_mass_and_deriv(z, Particles::tL)[0];
  const double mbsq = get_squared_mass_and_deriv(z, Particles::bL)[0];
  const double mhsq = get_squared_mass_and_deriv(z, Particles::h)[0];
  const double GSS  = 4.9e-4 * Ki->Tc;
  const double GY   = 4.2e-3 * Ki->Tc;
  const double GM   = mtsq / (63. * Ki->Tc);
  const double Gh =
      get_squared_mass_and_deriv(z, Particles::W)[0] / (50. * Ki->Tc);
  const double Gtott = K(K4, fermion, std::sqrt(mtsq)) * Ki->Tc /
                       (6. * K(D0, fermion, std::sqrt(mtsq)));
  const double Gtotb = K(K4, fermion, std::sqrt(mbsq)) * Ki->Tc /
                       (6. * K(D0, fermion, std::sqrt(mbsq)));
  const double Gtoth = K(K4, boson, std::sqrt(mhsq)) * Ki->Tc /
                       (20. * K(D0, fermion, std::sqrt(mhsq)));
  const double GW  = Gtoth;
  const double D0t = K(D0, fermion, std::sqrt(mtsq));
  const double D0b = K(D0, fermion, std::sqrt(mbsq));

  VecDoub Ci = {K0 * (GY + GW + GM + GSS * (1. + 9. * D0t)),
                0.,
                K0 * (-GW + GSS * (1. + 9. * D0b)),
                0.,
                -K0 * (GY + GM + GSS * (1. - 9. * D0t)),
                0.,
                K0 * GY,
                0.};
  res.push_back(Ci);
  Ci    = 0. * Ci;
  Ci[1] = -Gtott;
  res.push_back(Ci);
  Ci = {K0 * (-GW + GSS * (1. + 9. * D0t)),
        0.,
        K0 * (GY + GW + GSS * (1. + 9. * D0b)),
        0.,
        -K0 * (GY + GSS * (1. - 9. * D0t)),
        0.,
        K0 * GY,
        0.};
  res.push_back(Ci);
  Ci    = 0. * Ci;
  Ci[3] = -Gtotb;
  res.push_back(Ci);
  Ci = {-K0 * (GY + GM + GSS * (1. + 9. * D0t)),
        0.,
        -K0 * (GY + GSS * (1. + 9. * D0b)),
        0.,
        K0 * (2 * GY + GM + GSS * (1. - 9. * D0t)),
        0.,
        -2. * K0 * GY,
        0.};
  res.push_back(Ci);
  Ci    = 0. * Ci;
  Ci[5] = -Gtott;
  res.push_back(Ci);
  Ci = {K0 * GY, 0., K0 * GY, 0., -2. * K0 * GY, 0., K0 * (2. * GY + Gh), 0.};
  res.push_back(Ci);
  Ci    = 0. * Ci;
  Ci[7] = -Gtoth;
  res.push_back(Ci);
  return res;
}

VecDoub TransportNetwork::calc_Source(const double z)
{
  VecDoub res(2 * prtcl_list.size());
  const double dth  = calculate_theta(z, 1);
  const double d2th = calculate_theta(z, 2);
  VecDoub temp;
  double msq;
  double dmsq;
  Kfactor K(Ki, false);
  P_type pt;
  for (size_t i = 0; i < prtcl_list.size(); i++)
  {
    pt = get_particle_type(prtcl_list[i]);
    if ((pt == fermion) || (pt == antifermion))
    {
      temp       = get_squared_mass_and_deriv(z, prtcl_list[i]);
      msq        = temp[0];
      dmsq       = temp[1];
      double p1  = (dmsq * dth + msq * d2th);
      double p2  = dmsq * msq * dth;
      double hel = (pt == fermion) ? -1. : 1.;
      res[2 * i] =
          -Ki->vw * Ki->gamw * hel *
          (p1 * K(Q8o1, pt, std::sqrt(msq)) - p2 * K(Q9o1, pt, std::sqrt(msq)));
      res[2 * i + 1] =
          -Ki->vw * Ki->gamw * hel *
          (p1 * K(Q8o2, pt, std::sqrt(msq)) - p2 * K(Q9o2, pt, std::sqrt(msq)));
    }
    else
    {
      res[2 * i]     = 0.;
      res[2 * i + 1] = 0.;
    }
  }
  return res;
}

void TransportNetwork::operator()(const double z, VecDoub &u, VecDoub &du)
{
  /* MatDoub temp = calc_A_inv(z) * (calc_Collision(z));
  std::ofstream yfile("res.dat");
  for (auto it : temp)
  {
    for (auto jt : it)
      yfile << jt << "\t";
    yfile << "\n";
  }
  yfile.close();

  exit(1); */

  du = calc_A_inv(z) * ((calc_Collision(z) - calc_B(z)) * u + calc_Source(z));
}

} // namespace BSMPT
