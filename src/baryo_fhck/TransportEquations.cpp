#include <BSMPT/baryo_fhck/TransportEquations.h>

namespace BSMPT
{

namespace Baryo
{
namespace FHCK
{
TransportEquations::TransportEquations(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<CoexPhases> &CoexPhase_in,
    const double &vwall_in,
    const double &Tstar_in,
    const VevProfileMode &VevProfile_In)
{
  modelPointer = pointer_in;
  Tstar        = Tstar_in;
  vwall        = vwall_in;
  CoexPhase    = CoexPhase_in;
  FalseVacuum  = CoexPhase->false_phase.Get(Tstar).point;
  TrueVacuum   = CoexPhase->true_phase.Get(Tstar).point;
  VevProfile   = VevProfile_In;
  Initialize();
}

void TransportEquations::Initialize()
{
  EmptyVacuum = std::vector<double>(modelPointer->get_NHiggs(), 0);

  SetEtaInterface();

  if (VevProfile == VevProfileMode::FieldEquation)
  {
    const double &eps = 0.01;
    // Calculate tunnel profile
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      std::vector<double> res = modelPointer->MinimizeOrderVEV(vev);
      return modelPointer->VEff(res, Tstar);
    };

    std::function<std::vector<double>(std::vector<double>)> dV =
        [=](auto const &arg) { return NablaNumerical(arg, V, eps); };
    std::function<std::vector<std::vector<double>>(std::vector<double>)>
        Hessian = [=](auto const &arg)
    { return HessianNumerical(arg, V, eps); };

    vacuumprofile =
        std::make_unique<VacuumProfileNS::VacuumProfile>(FalseVacuum.size(),
                                                         TrueVacuum,
                                                         FalseVacuum,
                                                         V,
                                                         dV,
                                                         Hessian,
                                                         Lw.value());
    vacuumprofile->CalculateProfile();
  }

  zList = std::vector<double>(
      MakeDistribution(LwMultiplier * Lw.value(), NumberOfSteps));

  Logger::Write(LoggingLevel::FHCK, "Calculating fermion masses.\n");

  // Generate the fermion masses splines
  GenerateFermionMass();

  gamwall = 1. / std::sqrt(1. - vwall * vwall);

  nEqs = moment * (nFermions + nBosons);

  BuildKernelInterpolation();
}

tk::spline TransportEquations::InterpolateKernel(const std::string &kernel_name,
                                                 const bool is_1D)
{
  double x, x_save, vw, y;
  std::vector<double> vw_vals, y_vals, x_vals, res_vals;

  std::string filename = "kernels/" + kernel_name + ".dat";
  std::ifstream kernel(filename);
  if (is_1D)
  {
    while (kernel >> x >> y)
    {
      x_vals.push_back(x);
      res_vals.push_back(y);
    }
    kernel.close();
  }
  else
  {
    kernel >> x >> vw >> y;
    x_save = x;
    vw_vals.push_back(vw);
    y_vals.push_back(y);

    while (kernel >> x >> vw >> y)
    {
      if (x != x_save)
      {
        tk::spline spl(vw_vals, y_vals);
        x_vals.push_back(x_save);
        res_vals.push_back(spl(vwall));
        vw_vals.clear();
        y_vals.clear();
      }
      vw_vals.push_back(vw);
      y_vals.push_back(y);
      x_save = x;
    }
    kernel.close();
  }

  return tk::spline(x_vals, res_vals);
}

void TransportEquations::BuildKernelInterpolation()
{
  std::string kernel_name;
  bool is_1d = true;

  K4FHf = InterpolateKernel("K4FH_f", true);
  K4FHb = InterpolateKernel("K4FH_b", true);
  Rbarf = InterpolateKernel("Rbar_f", false);
  Rbarb = InterpolateKernel("Rbar_b", false);
  Qlf.push_back(tk::spline());
  Q8ol.push_back(tk::spline());
  Q9ol.push_back(tk::spline());
  Qlb.push_back(tk::spline());

  for (size_t l = 0; l <= moment; l++)
  {
    for (int type = 0; type <= 1; type++)
    {
      std::string suffix = (type == 0 ? "_f" : "_b");

      kernel_name = "D" + std::to_string(l) + suffix;
      type == 0 ? Dlf.push_back(InterpolateKernel(kernel_name, is_1d))
                : Dlb.push_back(InterpolateKernel(kernel_name, is_1d));

      if (l == 1)
        type == 0 ? Klf.push_back(tk::spline()) : Klb.push_back(tk::spline());
      else
      {
        kernel_name = "K" + std::to_string(l) + suffix;
        type == 0 ? Klf.push_back(InterpolateKernel(kernel_name, is_1d))
                  : Klb.push_back(InterpolateKernel(kernel_name, is_1d));
      }

      if (l != 0)
      {
        kernel_name = "Q" + std::to_string(l) + suffix;
        type == 0 ? Qlf.push_back(InterpolateKernel(kernel_name, is_1d))
                  : Qlb.push_back(InterpolateKernel(kernel_name, is_1d));

        if (type == 0)
        {
          kernel_name = "Q8o" + std::to_string(l) + suffix;
          Q8ol.push_back(InterpolateKernel(kernel_name, is_1d));

          kernel_name = "Q9o" + std::to_string(l) + suffix;
          Q9ol.push_back(InterpolateKernel(kernel_name, is_1d));
        }
      }
    }
    is_1d = false;
  }
}

void TransportEquations::SetNumberOfSteps(const int &num)
{
  NumberOfSteps = num;
  Initialize();
}

void TransportEquations::SetEtaInterface()
{
  auto config =
      std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
  EtaInterface =
      std::make_shared<CalculateEtaInterface>(config, GetSMConstants());

  EtaInterface->CalcEta(vwall,
                        TrueVacuum,
                        FalseVacuum,
                        Tstar,
                        modelPointer,
                        Minimizer::WhichMinimizerDefault);

  Lw = EtaInterface->getLW();

  Logger::Write(LoggingLevel::FHCK,
                "Lw * T = " + std::to_string(Lw.value() * Tstar) + "\n");
}

std::vector<double> TransportEquations::Vev(const double &z, const int &diff)
{
  if (VevProfile == VevProfileMode::Kink)
  {
    if (diff == 0)
      return FalseVacuum +
             (1 - tanh(z / Lw.value())) * (TrueVacuum - FalseVacuum) / 2;
    if (diff == 1)
      return -1 / pow(cosh(z / Lw.value()), 2) * (TrueVacuum - FalseVacuum) /
             (2 * Lw.value());
    if (diff == 2)
      return 2 * tanh(z / Lw.value()) / pow(cosh(z / Lw.value()), 2) *
             (TrueVacuum - FalseVacuum) / (2 * Lw.value() * Lw.value());
  }
  else if (VevProfile == VevProfileMode::FieldEquation)
  {
    return vacuumprofile->GetVev(z, diff);
  }

  Logger::Write(LoggingLevel::FHCK, "Error! No VEV profile selected. - Vev()");
  std::runtime_error("VEV profile mode selected is not valid.");
  return std::vector<double>();
}

void TransportEquations::GenerateFermionMass()
{
  QuarkMassesRe.clear();
  QuarkMassesIm.clear();
  size_t ind = (modelPointer->get_NQuarks() - 1);
  std::vector<std::vector<double>> MassesReal, MassesImag;
  for (auto z : zList)
  {
    std::vector<double> MassReal, MassImag;

    Eigen::MatrixXcd MIJQuarks =
        modelPointer->QuarkMassMatrix(modelPointer->MinimizeOrderVEV(Vev(z)));
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esQuark(MIJQuarks);

    for (auto m : esQuark.eigenvalues())
    {
      const int fac = 2 * (m.real() >= 0) - 1;
      MassReal.push_back(m.real() / fac);
      MassImag.push_back(m.imag() / fac);
    }
    MassesReal.push_back(MassReal);
    MassesImag.push_back(MassImag);
  }

  MassesReal = Transpose(MassesReal);
  MassesImag = Transpose(MassesImag);

  for (auto MassProfile : MassesReal)
    QuarkMassesRe.push_back(tk::spline(zList,
                                       MassProfile,
                                       tk::spline::cspline,
                                       false,
                                       tk::spline::not_a_knot,
                                       0,
                                       tk::spline::not_a_knot,
                                       0));

  for (auto MassProfile : MassesImag)
    QuarkMassesIm.push_back(tk::spline(zList,
                                       MassProfile,
                                       tk::spline::cspline,
                                       false,
                                       tk::spline::not_a_knot,
                                       0,
                                       tk::spline::not_a_knot,
                                       0));

  AsciiPlotter Plot("Top mass", 120, ceil(120 / 3.));

  Plot.addPlot(zList, MassesReal.at(ind), "Re(mt)", '*');
  Plot.addPlot(zList, MassesImag.at(ind), "Im(mt)", '.');
  Plot.legend();
  std::stringstream ss;
  Plot.show(ss);
  Logger::Write(LoggingLevel::VacuumProfile, ss.str());
}

void TransportEquations::GetFermionMass(const double &z,
                                        const size_t &fermion,
                                        double &m2,
                                        double &m2prime,
                                        double &thetaprime,
                                        double &theta2prime)
{
  std::complex<double> m, mprime, mprimeprime;
  const int ind =
      (modelPointer->get_NQuarks() - 1) - fermion; // Index of fermion

  m = std::complex<double>(QuarkMassesRe.at(ind)(z), QuarkMassesIm.at(ind)(z));

  mprime      = std::complex<double>(QuarkMassesRe.at(ind).deriv(1, z),
                                QuarkMassesIm.at(ind).deriv(1, z));
  mprimeprime = std::complex<double>(QuarkMassesRe.at(ind).deriv(2, z),
                                     QuarkMassesIm.at(ind).deriv(2, z));
  m2          = std::abs(m * m) / (Tstar * Tstar); // m^2 of fermion
  m2prime     = std::abs(2. * mprime * m) / (Tstar * Tstar);
  // Calculate theta
  if (VevProfile == VevProfileMode::Kink)
  {
    double brk;
    double sym;

    // Use BSMPTv2 functions
    switch (fermion)
    {
    case 0:
      brk = EtaInterface->getBrokenCPViolatingPhase_top();
      sym = EtaInterface->getSymmetricCPViolatingPhase_top();
      break;
    case 1:
      brk = EtaInterface->getBrokenCPViolatingPhase_bot();
      sym = EtaInterface->getSymmetricCPViolatingPhase_bot();
      break;
    case 2:
      brk = EtaInterface->getBrokenCPViolatingPhase_top();
      sym = EtaInterface->getSymmetricCPViolatingPhase_top();
      break;
    default: throw("Invalid fermion in GetFermionMass()"); break;
    }
    thetaprime =
        -0.5 * ((brk - sym) * 1 / pow(cosh(z / Lw.value()), 2)) / Lw.value();
    theta2prime =
        ((brk - sym) / pow(cosh(z / Lw.value()), 2) * tanh(z / Lw.value())) /
        pow(Lw.value(), 2);
  }
  else if (VevProfile == VevProfileMode::FieldEquation)
  {
    if (pow(pow(m.imag(), 2) + pow(m.real(), 2), 2) == 0) // avoid 1/0
    {
      thetaprime  = 0.;
      theta2prime = 0.;
    }
    else
    {
      // Deduce theta' and theta'' from the mass's phase
      thetaprime = (m.real() * mprime.imag() - m.imag() * mprime.real()) /
                   (pow(m.imag(), 2) + pow(m.real(), 2));
      theta2prime =
          (pow(m.real(), 2) * (-2 * mprime.imag() * mprime.real() +
                               m.real() * mprimeprime.imag()) +
           pow(m.imag(), 2) * (2 * mprime.imag() * mprime.real() +
                               m.real() * mprimeprime.imag()) -
           pow(m.imag(), 3) * mprimeprime.real() -
           m.imag() * m.real() *
               (2 * pow(mprime.imag(), 2) - 2 * pow(mprime.real(), 2) +
                m.real() * mprimeprime.real())) /
          pow(pow(m.imag(), 2) + pow(m.real(), 2), 2);
    }
  }
}

double TransportEquations::GetWMass(const std::vector<double> &vev,
                                    const double &T) const
{
  std::vector<double> res;
  res =
      modelPointer->GaugeMassesSquared(modelPointer->MinimizeOrderVEV(vev), T);
  std::vector<double> nrepeat(modelPointer->get_NGauge());
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
      return sqrt(res.at(j)) / Tstar;
    }
  }
  return 0;
}

MatDoub TransportEquations::CalculateCollisionMatrix(const double &mW,
                                                     VecDoub &FermionMasses,
                                                     VecDoub &BosonMasses)
{
  if (FermionMasses.size() != 3)
    throw std::runtime_error(
        "Transport equation are only supported with 3 Weyl fermions.");

  if (BosonMasses.size() != 1)
    throw std::runtime_error(
        "Transport equation are only supported with 1 boson.");

  std::vector<double> mass = {sqrt(FermionMasses[0]),
                              sqrt(FermionMasses[0]),
                              sqrt(FermionMasses[2]),
                              sqrt(BosonMasses[0])};

  std::vector<ParticleType> ptype = {Fermion, Fermion, Fermion, Boson};

  const double D0T = Dlf[0](mass[0]);
  const double D0B = Dlf[0](mass[2]);

  const std::vector<double> gammatot = {K4FHf(mass[0]) / (D0T * 6.) * Tstar,
                                        K4FHf(mass[1]) / (D0T * 6.) * Tstar,
                                        K4FHf(mass[2]) / (D0B * 6.) * Tstar,
                                        K4FHb(mass[3]) /
                                            (Dlb[0](mass[3]) * 20.) * Tstar};

  const std::vector<double> rates = {pow(mass[0], 2) / 63. * Tstar, // gammaM
                                     4.2e-3 * Tstar,                // gammaY
                                     gammatot[3],                   // gammaW
                                     mW * mW / 50. * Tstar,         // gammaH
                                     4.9e-4 * Tstar};               // gammaSS

  std::vector<double> trunc = {1., 1., 1., 1., 1.};

  const std::vector<std::vector<double>> tL_interactions = {
      {1., 1., 1., 0., (1. + 9. * D0T)},
      {-1., -1., 0., 0., (-1. + 9. * D0T)},
      {0., 0., -1., 0., (1. + 9. * D0B)},
      {0., 1., 0., 0., 0.}};
  const std::vector<std::vector<double>> tR_interactions = {
      {-1., -1., 0., 0., -(1. + 9. * D0T)},
      {1., 2., 0., 0., (1. - 9. * D0T)},
      {0., -1., 0., 0., -(1. + 9. * D0B)},
      {0., -2., 0., 0., 0.}};
  const std::vector<std::vector<double>> bL_interactions = {
      {0., 0., -1., 0., (1. + 9. * D0T)},
      {0., -1., 0., 0., (-1. + 9. * D0T)},
      {0., 1., 1., 0., (1. + 9. * D0B)},
      {0., 1., 0., 0., 0.}};
  const std::vector<std::vector<double>> h_interactions = {
      {0., 1., 0., 0., 0.},
      {0., -2., 0., 0., 0},
      {0., 1., 0., 0., 0.},
      {0., 2., 0., 1., 0.}};
  const std::vector<std::vector<std::vector<double>>> prtcl = {
      tL_interactions, tR_interactions, bL_interactions, h_interactions};

  MatDoub Gamma(nEqs, nEqs, 0.);
  double Kp, ul;
  for (size_t i = 0; i < prtcl.size(); i++)
  {
    double m = mass[i];
    for (size_t l = 1; l <= moment; l++)
    {
      if (l == 1)
      {
        ul                      = 0.;
        trunc[rates.size() - 1] = 1.;
        Kp = (ptype[i] == Fermion ? Klf[l - 1](m) : Klb[l - 1](m));
      }
      else if (l == 2)
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 1.;
        Kp = -vwall * (ptype[i] == Fermion ? Klf[0](m) : Klb[0](m));
      }
      else
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 0.;
        Kp = (ptype[i] == Fermion ? Klf[l - 1](m) : Klb[l - 1](m));
      }

      for (size_t j = 0; j < prtcl.size(); j++)
      {
        for (size_t k = 0; k < rates.size(); k++)
        {
          Gamma[moment * i + l - 1][moment * j] +=
              prtcl[i][j][k] * rates[k] * trunc[k];
        }
        Gamma[moment * i + l - 1][moment * j] *= Kp;
      }
      Gamma[moment * i + l - 1][moment * i + l - 1] -= ul * gammatot[i];
    }
  }

  return Gamma;
}

MatDoub TransportEquations::calc_Ainv(const double &m, const ParticleType &type)
{
  MatDoub res(moment, moment, 0.);
  std::vector<double> Di(moment, 1), Ri(moment);

  for (size_t i = 0; i < moment - 1; i++)
  {
    // 1, D1, D2, ..., Dn-1
    Di.at(i + 1) = (type == Fermion ? Dlf[i + 1](m) : Dlb[i + 1](m));
    // off-diagonal "diagonal"
    res[i + 1][i] = 1;
  }

  // 0, 0, ..., -vwall, -1
  Ri.at(moment - 2) = -vwall;
  Ri.at(moment - 1) = -1;

  // (-1)^moment Inv(A)
  double Dn = (type == Fermion ? Dlf[moment](m) : Dlb[moment](m));
  for (size_t i = 0; i < moment - 1; i++)
    Dn -= Ri.at(i) * Di.at(i + 1);

  // Tensor product
  for (std::size_t i = 0; i < moment; i++)
    for (std::size_t j = 0; j < moment; j++)
      res[i][j] += Di[i] * Ri[j] / Dn;

  return res;
}

MatDoub TransportEquations::calc_m2B(const double &m,
                                     const double &dm2,
                                     const ParticleType &type)
{
  MatDoub res(moment, moment, 0.);
  const double fRbar = (type == Fermion ? Rbarf(m) : Rbarb(m));

  for (size_t l = 1; l <= moment; l++)
  {
    const double fQ   = (type == Fermion ? Qlf[l](m) : Qlb[l](m));
    res[l - 1][l - 1] = fRbar * (double)(l - 1) * dm2;
    res[l - 1][0]     = gamwall * vwall * fQ * dm2;
  }
  return res;
}

VecDoub TransportEquations::calc_source(const double &m,
                                        const double &dm2,
                                        const double &dth,
                                        const double &d2th,
                                        const ParticleType &type)
{
  VecDoub res(moment);
  for (size_t i = 0; i < moment; i++)
    res[i] = -vwall * gamwall *
             ((dm2 * dth + m * m * d2th) * Q8ol[i + 1](m) -
              dm2 * m * m * dth * Q9ol[i + 1](m)); // Si
  return res;
}

void TransportEquations::Equations(const double &z,
                                   MatDoub &Mtilde,
                                   VecDoub &Stilde)
{
  const int nF2 = 2 * nFermions; // Clearer code

  MatDoub Ainverse(nEqs, nEqs, 0.); // Store A^-1

  VecDoub S(nEqs, 0.); // Source vector

  MatDoub m2B(nEqs, nEqs, 0.); // m2' B;

  Stilde.zero(); // A^-1 * Source vector
  Mtilde.zero(); // Store A^-1 * M

  // Quark matrix
  Eigen::MatrixXcd MIJQuarks =
      modelPointer->QuarkMassMatrix(modelPointer->MinimizeOrderVEV(Vev(z, 0)));
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esQuark(MIJQuarks);

  // Mass vector
  const double mW = GetWMass(Vev(z, 0), Tstar);
  VecDoub FermionMasses(nFermions);
  VecDoub BosonMasses(nBosons);

  // Fermions
  for (size_t fermion = 0; fermion < nFermions; fermion++)
  {

    double m2, m2prime, thetaprime, theta2prime;
    // Calculate fermionic masses
    GetFermionMass(z, fermion, m2, m2prime, thetaprime, theta2prime);
    FermionMasses[fermion] = m2;

    // Calculate A inverse for fermion
    MatDoub tempA(calc_Ainv(sqrt(m2), ParticleType::Fermion));
    Ainverse[2 * fermion][2 * fermion]         = tempA[0][0];
    Ainverse[2 * fermion][2 * fermion + 1]     = tempA[0][1];
    Ainverse[2 * fermion + 1][2 * fermion]     = tempA[1][0];
    Ainverse[2 * fermion + 1][2 * fermion + 1] = tempA[1][1];

    // Calculate m2'B for fermion
    MatDoub tempB(calc_m2B(sqrt(m2), m2prime, ParticleType::Fermion));
    m2B[2 * fermion][2 * fermion]         = tempB[0][0];
    m2B[2 * fermion + 1][2 * fermion]     = tempB[1][0];
    m2B[2 * fermion + 1][2 * fermion + 1] = tempB[1][1];

    // Source terms
    VecDoub tempS(calc_source(
        sqrt(m2), m2prime, thetaprime, theta2prime, ParticleType::Fermion));
    S[2 * fermion]     = tempS[0];
    S[2 * fermion + 1] = tempS[1];
  }

  // Bosons
  for (size_t boson = 0; boson < nBosons; boson++)
  {
    double m2, m2prime;
    // Calculate the boson mass
    m2                 = 0;
    m2prime            = 0;
    BosonMasses[boson] = m2;

    // Calculate A inverse for bosons
    MatDoub tempA(calc_Ainv(sqrt(m2), ParticleType::Boson));
    Ainverse[nF2 + 2 * boson][nF2 + 2 * boson]         = tempA[0][0];
    Ainverse[nF2 + 2 * boson][nF2 + 2 * boson + 1]     = tempA[0][1];
    Ainverse[nF2 + 2 * boson + 1][nF2 + 2 * boson]     = tempA[1][0];
    Ainverse[nF2 + 2 * boson + 1][nF2 + 2 * boson + 1] = tempA[1][1];

    // Calculate m2'B
    MatDoub tempB(calc_m2B(sqrt(m2), m2prime, ParticleType::Fermion));
    m2B[nF2 + 2 * boson][nF2 + 2 * boson]         = tempB[0][0];
    m2B[nF2 + 2 * boson + 1][nF2 + 2 * boson]     = tempB[1][0];
    m2B[nF2 + 2 * boson + 1][nF2 + 2 * boson + 1] = tempB[1][1];
  }

  // Gamma = deltaC - m2' B
  const MatDoub CollisiontMatrix =
      CalculateCollisionMatrix(mW, FermionMasses, BosonMasses);

  // Calculate M = A^-1 * Gamma ( = deltaC - m2'B)
  for (size_t i = 0; i < 2 * (nBosons + nFermions); i++)
    for (size_t j = 0; j < 2 * (nBosons + nFermions); j++)
      for (size_t l = 0; l < 2 * (nBosons + nFermions); l++)
        Mtilde[i][j] += Ainverse[i][l] * (CollisiontMatrix[l][j] - m2B[l][j]);

  // Calculate Stilde = A^-1 * S
  for (size_t i = 0; i < 2 * (nBosons + nFermions); i++)
    for (size_t j = 0; j < 2 * (nFermions); j++)
      Stilde[i] += Ainverse[i][j] * S[j];
}

std::vector<double> TransportEquations::MakeDistribution(const double xmax,
                                                         const size_t npoints)
{
  std::vector<double> res(npoints);

  for (size_t i = 0; i < npoints; i++)
  {
    double temp = pow(((i - npoints / 2.)) / (npoints / 2.), 3) * M_PI / 4.;
    res.at(i)   = xmax * tan(temp);
  }
  return res;
}

void TransportEquations::CheckBoundary(const MatDoub &MtildeM,
                                       const VecDoub &StildeM,
                                       const MatDoub &MtildeP,
                                       const VecDoub &StildeP)
{
  double STildeLength(0);
  size_t NumberOfNonDecayingModes(0);
  stringstream ss;

  Eigen::MatrixXcd EigenMtildeM(nEqs, nEqs);
  Eigen::MatrixXcd EigenMtildeP(nEqs, nEqs);

  for (size_t i = 0; i < nEqs; i++)
  {
    STildeLength += pow(StildeM[i], 2);
    STildeLength += pow(StildeP[i], 2);
    for (size_t j = 0; j < nEqs; j++)
    {
      EigenMtildeM(i, j) = MtildeM[i][j];
      EigenMtildeP(i, j) = MtildeP[i][j];
    }
  }

  // Checking if S vector are small enough
  STildeLength = sqrt(STildeLength);
  STildeLength /= 2. * nEqs;

  if (STildeLength > STildeThreshold)
  {
    Status = FHCKStatus::SmallIntegrationRegion;
    ss << "\033[31mThe integration region is too small, increase "
          "LwMultiplier.\033[0m\n";
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> EigenSolverM(EigenMtildeM);
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> EigenSolverP(EigenMtildeP);

  // Number of modes that have to be set to zero
  for (auto ev : EigenSolverM.eigenvalues())
    if (ev.real() >= 0) NumberOfNonDecayingModes++;
  for (auto ev : EigenSolverP.eigenvalues())
    if (ev.real() >= 0) NumberOfNonDecayingModes++;

  ss << "There are " << NumberOfNonDecayingModes
     << " modes that have to be set to zero at the boundaries.";

  if (NumberOfNonDecayingModes > nEqs)
  {
    ss << " \033[31m\nToo many non-decaying modes. Impossible to satisfy "
          "the "
          "boundary "
          "conditions of mu = u = 0.\033[0m\n";
    Status = FHCKStatus::UnphysicalBoundary;
  }
  else
    ss << " \033[92mSuccess.\033[0m";

  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void TransportEquations::SolveTransportEquation()
{
  Logger::Write(LoggingLevel::FHCK, "Lw = " + std::to_string(Lw.value()));

  MatDoub STildeList(NumberOfSteps, nEqs);
  Mat3DDoub MTildeList(NumberOfSteps, nEqs, nEqs);

  MatDoub Mtilde(nEqs, nEqs), MtildeM(nEqs, nEqs), MtildeP(nEqs, nEqs);
  VecDoub Stilde(nEqs), StildeM(nEqs), StildeP(nEqs);

  Equations(-1, MtildeM, StildeM);
  Equations(1, MtildeP, StildeP);

  CheckBoundary(MtildeM, StildeM, MtildeP, StildeP);
  if (Status != FHCKStatus::NotSet)
  {
    Logger::Write(
        LoggingLevel::FHCK,
        "\033[31mTransport equation unable to be computed. Status code : " +
            FHCKStatusToString.at(FHCKStatus::NotSet) + "\033[0m\n");
    return;
  }

  for (size_t i = 1; i < NumberOfSteps; i++)
  {
    // Compute Mtilde and Stilde
    double zc = (zList[i] + zList[i - 1]) / 2.;
    if (zc < -1) // TODO: Fix this, too specific.
    {
      Mtilde = MtildeM;
      Stilde = StildeM;
    }
    else if (zc > 1)
    {
      Mtilde = MtildeP;
      Stilde = StildeP;
    }
    else
      Equations(zc, Mtilde, Stilde);

    // Save the Mtilde and Stilde
    for (size_t j = 0; j < nEqs; j++)
    {
      STildeList[i][j] = Stilde[j];
      for (size_t k = 0; k < nEqs; k++)
        MTildeList[i][j][k] = Mtilde[j][k];
    }
  }
  // Construct Difeq object (S_j,n matrix)
  Difeq_TransportEquation difeq(
      zList, nFermions, nBosons, MTildeList, STildeList, moment);

  size_t itmax = 30;
  double conv  = 1e-10;
  double slowc = 1e-3;
  VecDoub scalv(zList.size(), 1);
  VecInt indexv(nEqs);
  indexv[0] = 0;
  indexv[1] = 4;
  indexv[2] = 1;
  indexv[3] = 5;
  indexv[4] = 2;
  indexv[5] = 6;
  indexv[6] = 3;
  indexv[7] = 7;
  int NB    = 4;
  MatDoub y(nEqs, zList.size(), 0.);
  RelaxOde solvde(itmax, conv, slowc, scalv, indexv, NB, y, difeq);

  std::string str = "test.csv";
  std::ofstream res(str);
  for (size_t i = 0; i < zList.size(); i++)
  {
    res << zList[i] << "\t" << y[0][i] << "\t" << y[1][i] << "\t" << y[2][i]
        << "\t" << y[3][i] << "\t" << y[4][i] << "\t" << y[5][i] << "\t"
        << y[6][i] << "\t" << y[7][i] << "\n";
  }
  res.close();

  // Store the solution
  Solution = y;

  PrintTransportEquation(120, "tL", "mu", 100);
  PrintTransportEquation(120, "tR", "mu", 100);
  PrintTransportEquation(120, "bL", "mu", 100);
  PrintTransportEquation(120, "h", "mu", 100);

  CalculateBAU();
}

void TransportEquations::CalculateBAU()
{
  double mt2, mb2, m2prime, thetaprime, theta2prime; // temporary vars
  // Weak spharelon rate
  const double Gsph = 1.e-6 * Tstar;
  double r;                // temporary variable to store the result
  std::vector<double> z;   // list of z positions
  std::vector<double> muB; // muB integrand at position z
  for (size_t i = 0; i < zList.size(); i++)
  {

    const double zi = zList[i]; // z at position i
    // Calculate the vev at z
    const std::vector<double> vev = Vev(zi);
    GetFermionMass(zi, 0, mt2, m2prime, thetaprime, theta2prime);
    GetFermionMass(zi, 2, mb2, m2prime, thetaprime, theta2prime);
    // Results
    r = 0;
    r += (1 + 4 * Dlf[0](sqrt(mt2))) / 2. * Solution.value()[0][i]; // tL
    r += (1 + 4 * Dlf[0](sqrt(mb2))) / 2. * Solution.value()[4][i]; // bL
    r += 2. * Dlf[0](sqrt(mt2)) * Solution.value()[2][i];           // tR
    r *= min(
        1.,
        2.4 * Tstar / Gsph *
            exp(-40 *
                modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(vev), 0.) /
                Tstar)); // f_sph(z)
    r *= exp(-45 * Gsph * std::abs(zi) /
             (4. * vwall * gamwall)); // exp(-45 G_sph |z| / 4 vw gammaw)
    // Save in list to pass to integrator
    z.push_back(zi);
    muB.push_back(r); // integrand
  }
  // Step 1: Set up GSL interpolation
  size_t n = z.size();
  const gsl_interp_type *interp_type =
      gsl_interp_cspline; // Cubic spline interpolation
  gsl_interp *interp = gsl_interp_alloc(interp_type, n);
  gsl_interp_init(interp, z.data(), muB.data(), n);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  // Step 2: Integrate the interpolated function
  double a = z.front(); // Lower bound of integration
  double b = z.back();  // Upper bound of integration
  double result, error;
  // Struct to hold parameters
  struct Params
  {
    gsl_interp *interp;
    gsl_interp_accel *acc;
    const std::vector<double> *z;
    const std::vector<double> *muB;
  } params                             = {interp, acc, &z, &muB};
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = [](double x, void *p) -> double
  {
    auto *data = static_cast<Params *>(p);
    return gsl_interp_eval(
        data->interp, data->z->data(), data->muB->data(), x, data->acc);
  };
  F.params = &params;
  gsl_integration_qags(
      &F, a, b, 1e-30, 1e-30, 1000, workspace, &result, &error);
  // Cleanup
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  // Results
  const double prefactor =
      405. * Gsph /
      (4 * M_PI * M_PI * vwall * gamwall * 106.75 * Tstar); // TODO Fix gstas

  result *= prefactor;
  // Save the result
  BAUEta = result;
  // Step 3: Output the result
  stringstream ss;
  ss << "eta = " << result << " with error " << error << std::endl;
  ss << "eta/eta_obs = " << result / (6.2e-10) << std::endl;
  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void TransportEquations::PrintTransportEquation(const int &size,
                                                const std::string &Particle,
                                                const std::string &MuOrU,
                                                const double &multiplier)
{
  AsciiPlotter Plot(Particle + " " + MuOrU, size, ceil(size / 3.));
  std::optional<int> ind;
  std::vector<double> z, y;

  if (Particle == "tL") ind = 0;
  if (Particle == "tR") ind = 2;
  if (Particle == "bL") ind = 4;
  if (Particle == "h") ind = 6;

  if (not ind.has_value()) throw("Invalid particle to plot the solution.");

  if (MuOrU == "u") ind = ind.value() + 1;

  for (size_t i = 0; i < zList.size(); i++)
    if (abs(zList[i]) < Lw.value() * multiplier)
    {
      z.push_back(zList[i]);
      y.push_back(Solution.value()[ind.value()][i]);
    }
  Plot.addPlot(z, y, "", '*');
  std::stringstream ss;
  Plot.show(ss);
  Logger::Write(LoggingLevel::FHCK, ss.str());
}
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT