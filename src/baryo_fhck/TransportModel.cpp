#include <BSMPT/baryo_fhck/TransportModel.h>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

TransportModel::TransportModel(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::vector<double> TrueVacuum_In,
    const std::vector<double> FalseVacuum_In,
    const double &vwall_in,
    const double &Tstar_in,
    const VevProfileMode &VevProfile_in,
    const TruncationScheme &truncationscheme_in)
{
  modelPointer     = pointer_in;
  Tstar            = Tstar_in;
  vwall            = vwall_in;
  TrueVacuum       = TrueVacuum_In;
  FalseVacuum      = FalseVacuum_In;
  VevProfile       = VevProfile_in;
  truncationscheme = truncationscheme_in;

  stringstream ss;
  ss << "----------------- Baryon Asymmetry Calculation -----------------\n";
  ss << "TrueVacuum = " << TrueVacuum << "\n";
  ss << "FalseVacuum = " << FalseVacuum << "\n";
  ss << "T = " << Tstar << "\n";
  ss << "vw = " << vwall << "\n";
  ss << "VEV profile = " << VevProfileModeToString.at(VevProfile) << "\n";
  ss << "Truncation scheme = " << TruncationSchemeToString.at(truncationscheme);

  Logger::Write(LoggingLevel::FHCK, ss.str());
}

TransportModel::TransportModel(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<CoexPhases> &CoexPhase,
    const double &vwall_in,
    const double &Tstar_in,
    const VevProfileMode &VevProfile_in,
    const TruncationScheme &truncationscheme_in)
    : TransportModel(pointer_in,
                     CoexPhase->true_phase.Get(Tstar_in).point,
                     CoexPhase->false_phase.Get(Tstar_in).point,
                     vwall_in,
                     Tstar_in,
                     VevProfile_in,
                     truncationscheme_in)
{
}

void TransportModel::Initialize()
{
  EmptyVacuum = std::vector<double>(modelPointer->get_NHiggs(), 0);

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

    vacuumprofile = std::make_unique<VacuumProfileNS::VacuumProfile>(
        FalseVacuum.size(), TrueVacuum, FalseVacuum, V, dV, Hessian);
    vacuumprofile->CalculateProfile();
    if (vacuumprofile->status != VacuumProfileNS::VacuumProfileStatus::Success)
    {
      Logger::Write(LoggingLevel::FHCK, "Vacuum Profile Calculation failed!");
      status = TransportModelStatus::Failed;
      return;
    }

    status = TransportModelStatus::Success;
    Lw     = vacuumprofile->Lw;
  }

  if (VevProfile == VevProfileMode::Kink)
  {
    SetEtaInterface();
  }
  status = TransportModelStatus::Success;
}

void TransportModel::SetEtaInterface()
{
  auto config =
      std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
  EtaInterface =
      std::make_shared<CalculateEtaInterface>(config, GetSMConstants());

  Logger::Write(LoggingLevel::FHCK, "Calculating Lw (EtaInterface)...");

  try
  {
    EtaInterface->CalcEta(vwall,
                          TrueVacuum,
                          FalseVacuum,
                          Tstar,
                          modelPointer,
                          Minimizer::WhichMinimizerDefault);
  }
  catch (...)
  {
    Logger::Write(LoggingLevel::FHCK, "Failed.");
    status = TransportModelStatus::Failed;
    return;
  }

  Logger::Write(LoggingLevel::FHCK, "\033[92mSuccess.\033[0m");
  status = TransportModelStatus::Success;
  Lw     = EtaInterface->getLW();

  Logger::Write(LoggingLevel::FHCK,
                "Lw * T = " + std::to_string(Lw * Tstar) + "\n");
}

void TransportModel::GenerateFermionMass(const std::vector<double> &zList)
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

  size_t i_left = 0, i_right = zList.size() - 1;

  double max_dmdz = -1;

  // calculate maximum
  for (const auto &zi : zList)
  {
    const double dmdz = sqrt(pow(QuarkMassesRe.at(ind).deriv(1, zi), 2) +
                             pow(QuarkMassesIm.at(ind).deriv(1, zi), 2));
    max_dmdz          = std::max(dmdz, max_dmdz);
  }

  for (size_t i = i_left; i <= i_right; i++)
  {
    const double zi   = zList[i];
    const double dmdz = sqrt(pow(QuarkMassesRe.at(ind).deriv(1, zi), 2) +
                             pow(QuarkMassesIm.at(ind).deriv(1, zi), 2));
    if (dmdz > max_dmdz / 100.)
    {
      i_left = i;
      break;
    }
  }

  for (size_t i = i_right; i >= i_left; i--)
  {
    const double zi   = zList[i];
    const double dmdz = sqrt(pow(QuarkMassesRe.at(ind).deriv(1, zi), 2) +
                             pow(QuarkMassesIm.at(ind).deriv(1, zi), 2));
    if (dmdz > max_dmdz / 100.)
    {
      i_right = i;
      break;
    }
  }

  std::vector<double> zListPlot, RePlot, ImPlot;
  for (size_t i = i_left; i <= i_right; i++)
  {
    zListPlot.push_back(zList[i]);
    RePlot.push_back(MassesReal.at(ind)[i]);
    ImPlot.push_back(MassesImag.at(ind)[i]);
  }

  Plot.addPlot(zListPlot, RePlot, "Re(mt)", '*');
  Plot.addPlot(zListPlot, ImPlot, "Im(mt)", '.');

  Plot.legend();
  std::stringstream ss;
  Plot.show(ss);
  Logger::Write(LoggingLevel::VacuumProfile, ss.str());
}

std::vector<double> TransportModel::Vev(const double &z, const int &diff)
{
  if (VevProfile == VevProfileMode::Kink)
  {
    if (diff == 0)
      return FalseVacuum + (1 - tanh(z / Lw)) * (TrueVacuum - FalseVacuum) / 2;
    if (diff == 1)
      return -1 / pow(cosh(z / Lw), 2) * (TrueVacuum - FalseVacuum) / (2 * Lw);
    if (diff == 2)
      return 2 * tanh(z / Lw) / pow(cosh(z / Lw), 2) *
             (TrueVacuum - FalseVacuum) / (2 * Lw * Lw);
  }
  else if (VevProfile == VevProfileMode::FieldEquation)
  {
    return vacuumprofile->GetVev(z, diff);
  }

  Logger::Write(LoggingLevel::FHCK, "Error! No VEV profile selected. - Vev()");
  std::runtime_error("VEV profile mode selected is not valid.");
  return std::vector<double>();
}

void TransportModel::GetFermionMass(const double &z,
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
      brk = EtaInterface->getBrokenCPViolatingPhase_top();
      sym = EtaInterface->getSymmetricCPViolatingPhase_top();
      break;
    case 2:
      brk = EtaInterface->getBrokenCPViolatingPhase_bot();
      sym = EtaInterface->getSymmetricCPViolatingPhase_bot();
      break;
    default: throw("Invalid fermion in GetFermionMass()"); break;
    }
    thetaprime = -0.5 * ((brk - sym) * 1 / pow(cosh(z / Lw), 2)) / Lw;
    theta2prime =
        ((brk - sym) / pow(cosh(z / Lw), 2) * tanh(z / Lw)) / pow(Lw, 2);
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

double TransportModel::GetWMass(const double &z, const double &T)
{
  std::vector<double> res;
  const std::vector<double> vev = Vev(z);
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

double TransportModel::EWSBVEV(const double &z)
{
  const std::vector<double> vev = Vev(z);
  return modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(vev), 0.);
}

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT