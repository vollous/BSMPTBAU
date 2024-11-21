#include <BSMPT/baryo_fhck/TransportEquations.h>

namespace BSMPT
{

namespace Baryo
{
namespace FHCK
{
TransportEquations::TransportEquations(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const double &vwall_in,
    const double &Tstar_in,
    const std::vector<double> &FalseVacuum_in,
    const std::vector<double> &TrueVacuum_in)
{
  modelPointer = pointer_in;
  Tstar        = Tstar_in;
  vwall        = vwall_in;
  FalseVacuum  = FalseVacuum_in;
  TrueVacuum   = TrueVacuum_in;
  VevProfile   = VevProfileMode::Kink;
  Initialize();
}

TransportEquations::TransportEquations(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<CoexPhases> &CoexPhase_in,
    const double &vwall_in,
    const double &Tstar_in)
{
  modelPointer = pointer_in;
  Tstar        = Tstar_in;
  vwall        = vwall_in;
  CoexPhase    = CoexPhase_in;
  FalseVacuum  = CoexPhase->false_phase.Get(Tstar).point;
  TrueVacuum   = CoexPhase->true_phase.Get(Tstar).point;
  VevProfile   = VevProfileMode::Kink;
  Initialize();
}

TransportEquations::TransportEquations(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const double &vwall_in,
    const double &Tstar_in,
    const std::shared_ptr<BounceActionInt> &ActionInt_in,
    const VevProfileMode &Mode_in)
{
  modelPointer = pointer_in;
  Tstar        = Tstar_in;
  vwall        = vwall_in;
  ActionInt    = ActionInt_in;
  VevProfile   = Mode_in;
  Initialize();
}

TransportEquations::TransportEquations(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<BounceSolution> &Bounce_in,
    const int &which_transition_temp,
    const VevProfileMode &Mode_in)
{
  modelPointer = pointer_in;
  Bounce       = Bounce_in;
  Tstar        = Bounce->CalcTransitionTemp(which_transition_temp);
  vwall        = Bounce->vwall;
  VevProfile   = Mode_in;
  Initialize();
}

TransportEquations::TransportEquations(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<GravitationalWave> &GW_in,
    const VevProfileMode &Mode_in)
{
  modelPointer = pointer_in;
  GW           = GW_in;
  VevProfile   = Mode_in;
  Initialize();
}

void TransportEquations::Initialize()
{
  EmptyVacuum = std::vector<double>(modelPointer->get_NHiggs(), 0);

  if (VevProfile == VevProfileMode::Kink) SetEtaInterface();

  Ki   = std::make_unique<Kinfo>(Tstar, vwall);
  Kfac = std::make_unique<Kfactor>(Ki, true);
}

void TransportEquations::SetNumberOfSteps(const int &num)
{
  NumberOfSteps = num;
  Initialize();
}

void TransportEquations::SetEtaInterface()
{
  if (VevProfile != VevProfileMode::Kink)
    throw("Cannot set theta is a non kink vev profile");

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
  else if (VevProfile == VevProfileMode::TunnelPath)
  {
    return std::vector<double>(2, 0.0);
  }

  Logger::Write(LoggingLevel::FHCK, "Error! No VEV profile selected. - Vev()");
  std::runtime_error("VEV profile mode selected is not valid.");
  return std::vector<double>();
}

void TransportEquations::GetFermionMass(
    const double &z,
    const int &fermion,
    const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> &esQuark,
    double &m2,
    double &m2prime,
    double &thetaprime,
    double &theta2prime)
{
  std::complex<double> m, mprime, mprimeprime;
  int ind = (modelPointer->get_NQuarks() - 1) - fermion; // Index of fermion
  m       = esQuark.eigenvalues()[ind];
  mprime  = esQuark.eigenvectors().col(ind).dot(
               (modelPointer->QuarkMassMatrix(
                    modelPointer->MinimizeOrderVEV(Vev(z, 1))) -
                modelPointer->QuarkMassMatrix(
                    modelPointer->MinimizeOrderVEV(EmptyVacuum))) *
               esQuark.eigenvectors().col(ind)) /
           esQuark.eigenvectors().col(ind).dot(esQuark.eigenvectors().col(ind));
  mprimeprime =
      esQuark.eigenvectors().col(ind).dot(
          (modelPointer->QuarkMassMatrix(
               modelPointer->MinimizeOrderVEV(Vev(z, 1))) -
           modelPointer->QuarkMassMatrix(
               modelPointer->MinimizeOrderVEV(EmptyVacuum))) *
          esQuark.eigenvectors().col(ind)) /
      esQuark.eigenvectors().col(ind).dot(esQuark.eigenvectors().col(ind));
  m2      = std::abs(m * m); // m^2 of fermion
  m2prime = std::abs(2. * mprime * m);
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
  else if (VevProfile == VevProfileMode::TunnelPath)
  {
    // Deduce theta' and theta'' from the mass's phase
    thetaprime = (m.real() * mprime.imag() - m.imag() * mprime.real()) /
                 (pow(m.imag(), 2) + pow(m.real(), 2));
    theta2prime = (pow(m.real(), 2) * (-2 * mprime.imag() * mprime.real() +
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
      return sqrt(res.at(j));
    }
  }
  return 0;
}

std::vector<std::vector<double>> TransportEquations::CalculateCollisionMatrix(
    const double &mW,
    const std::vector<double> &FermionMasses,
    const std::vector<double> &BosonMasses)
{
  if (FermionMasses.size() != 3)
    throw std::runtime_error(
        "Transport equation are only supported with 3 Weyl fermions.");

  if (BosonMasses.size() != 1)
    throw std::runtime_error(
        "Transport equation are only supported with 1 boson.");

  const double mTop   = sqrt(FermionMasses.at(0));
  const double mBot   = sqrt(FermionMasses.at(2));
  const double mHiggs = sqrt(BosonMasses.at(0));

  const double D0t = (*Kfac)(K_type::D0, P_type::fermion, mTop);
  const double D0b = (*Kfac)(K_type::D0, P_type::fermion, mBot);

  const double K0t = Tstar; // TODO: fix
  const double K0b = Tstar; // TODO: fix
  const double K0h = Tstar; // TODO: fix

  const double GammaTotTop =
      (*Kfac)(K_type::K4, P_type::fermion, mTop) /
      ((*Kfac)(K_type::D0, P_type::fermion, mTop) * (6 / Tstar)); // TODO: fix
  const double GammaTotBot =
      (*Kfac)(K_type::K4, P_type::fermion, mBot) /
      ((*Kfac)(K_type::D0, P_type::fermion, mBot) * (6 / Tstar)); // TODO: fix
  const double GammaTotHiggs =
      (*Kfac)(K_type::K4, P_type::boson, mHiggs) /
      ((*Kfac)(K_type::D0, P_type::boson, mHiggs) * (20 / Tstar)); // TODO: fix

  const double GammaM      = pow(mTop, 2) / (63 * Tstar);
  const double GammaY      = 4.2e-3 * Tstar;
  const double GammaW      = GammaTotHiggs;
  const double GammaTildeY = 4.9e-4 * Tstar;
  const double GammaH      = mW * mW / (50 * Tstar); // TODO: fix

  std::vector<std::vector<double>> Gamma = {
      {GammaM + GammaW + GammaY + (1 + 9 * D0t) * GammaY,
       0,
       -GammaM - GammaY + (-1 + 9 * D0t) * GammaY,
       0,
       -GammaW + (1 + 9 * D0b) * GammaY,
       0,
       GammaY,
       0},
      {-((GammaM + GammaW + GammaY + (1 + 9 * D0t) * GammaY) * K0t * vwall),
       -GammaTotTop,
       -((-GammaM - GammaY + (-1 + 9 * D0t) * GammaY) * K0t * vwall),
       0,
       -((-GammaW + (1 + 9 * D0b) * GammaY) * K0t * vwall),
       0,
       -(GammaY * K0t * vwall),
       0},
      {-GammaM - GammaY - (1 + 9 * D0t) * GammaY,
       0,
       GammaM + 2 * GammaY - (-1 + 9 * D0t) * GammaY,
       0,
       -GammaY - (1 + 9 * D0b) * GammaY,
       0,
       -2 * GammaY,
       0},
      {-((-GammaM - GammaY - (1 + 9 * D0t) * GammaY) * K0t * vwall),
       0,
       -((GammaM + 2 * GammaY - (-1 + 9 * D0t) * GammaY) * K0t * vwall),
       -GammaTotTop,
       -((-GammaY - (1 + 9 * D0b) * GammaY) * K0t * vwall),
       0,
       2 * GammaY * K0t * vwall,
       0},
      {-GammaW + (1 + 9 * D0t) * GammaY,
       0,
       -GammaY + (-1 + 9 * D0t) * GammaY,
       0,
       GammaW + GammaY + (1 + 9 * D0b) * GammaY,
       0,
       GammaY,
       0},
      {-((-GammaW + (1 + 9 * D0t) * GammaY) * K0b * vwall),
       0,
       -((-GammaY + (-1 + 9 * D0t) * GammaY) * K0b * vwall),
       0,
       -((GammaW + GammaY + (1 + 9 * D0b) * GammaY) * K0b * vwall),
       -GammaTotBot,
       -(GammaY * K0b * vwall),
       0},
      {GammaTildeY,
       0,
       -2 * GammaTildeY,
       0,
       GammaTildeY,
       0,
       GammaH + 2 * GammaTildeY,
       0},
      {-(GammaTildeY * K0h * vwall),
       0,
       2 * GammaTildeY * K0h * vwall,
       0,
       -(GammaTildeY * K0h * vwall),
       0,
       -((GammaH + 2 * GammaTildeY) * K0h * vwall),
       -GammaTotHiggs}};
  return Gamma;
}

void TransportEquations::Equations(const double &z,
                                   std::vector<std::vector<double>> &Mtilde,
                                   std::vector<double> &Stilde)
{
  const int nF2  = 2 * nFermions;             // Clearer code
  const int nFB2 = 2 * (nBosons + nFermions); // Twice number of particles

  std::vector<std::vector<double>> Ainverse(nFB2,
                                            vector<double>(nFB2)); // Store A^-1

  std::vector<double> S(nFB2, 0); // Source vector

  std::vector<std::vector<double>> m2B(nFB2, vector<double>(nFB2, 0)); // m2' B;

  Stilde = std::vector<double>(nFB2, 0); // A^-1 * Source vector
  Mtilde =
      std::vector<std::vector<double>>(nFB2,
                                       vector<double>(nFB2)); // Store A^1 * M
  // Quark matrix
  Eigen::MatrixXcd MIJQuarks =
      modelPointer->QuarkMassMatrix(modelPointer->MinimizeOrderVEV(Vev(z, 0)));
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esQuark(MIJQuarks);

  // Mass vector
  const double mW = GetWMass(Vev(z, 0), Tstar);
  std::vector<double> FermionMasses;
  std::vector<double> BosonMasses;

  // Fermions
  for (int fermion = 0; fermion < nFermions; fermion++)
  {

    double m2, m2prime, thetaprime, theta2prime;
    // Calculate fermionic masses
    GetFermionMass(z, fermion, esQuark, m2, m2prime, thetaprime, theta2prime);
    FermionMasses.push_back(m2);
    // Set variables
    const double fD1      = (*Kfac)(K_type::D1, P_type::fermion, sqrt(m2));
    const double fD2      = (*Kfac)(K_type::D2, P_type::fermion, sqrt(m2));
    const double fR       = -Ki->vw;
    const double fRbar    = (*Kfac)(K_type::Rbar, P_type::fermion, sqrt(m2));
    const double fQ1      = (*Kfac)(K_type::Q1, P_type::fermion, sqrt(m2));
    const double fQ2      = (*Kfac)(K_type::Q2, P_type::fermion, sqrt(m2));
    const double fm2prime = m2prime;
    const double fQ8o1    = (*Kfac)(K_type::Q8o1, P_type::fermion, sqrt(m2));
    const double fQ8o2    = (*Kfac)(K_type::Q8o2, P_type::fermion, sqrt(m2));
    const double fQ9o1    = (*Kfac)(K_type::Q9o1, P_type::fermion, sqrt(m2));
    const double fQ9o2    = (*Kfac)(K_type::Q9o2, P_type::fermion, sqrt(m2));

    // Calculate A inverse for fermion
    Ainverse[2 * fermion][2 * fermion]         = fR / (fD2 - fD1 * fR);
    Ainverse[2 * fermion][2 * fermion + 1]     = -1 / (fD2 - fD1 * fR);
    Ainverse[2 * fermion + 1][2 * fermion]     = fD2 / (fD2 - fD1 * fR);
    Ainverse[2 * fermion + 1][2 * fermion + 1] = -fD1 / (fD2 - fD1 * fR);
    // Source terms
    S[2 * fermion] = Ki->vw * Ki->gamw *
                     ((m2prime * thetaprime + m2 * theta2prime) * fQ8o1 -
                      m2prime * m2 * thetaprime * fQ9o1); // S1

    S[2 * fermion + 1] = Ki->vw * Ki->gamw *
                         ((m2prime * thetaprime + m2 * theta2prime) * fQ8o2 -
                          m2prime * m2 * thetaprime * fQ9o2); // S2
    // Calculate Gamma = deltaC - m2'B
    m2B[2 * fermion][2 * fermion]         = Ki->gamw * Ki->vw * fQ1 * fm2prime;
    m2B[2 * fermion + 1][2 * fermion]     = Ki->gamw * Ki->vw * fQ2 * fm2prime;
    m2B[2 * fermion + 1][2 * fermion + 1] = fRbar * fm2prime;
  }

  // Bosons
  for (int boson = 0; boson < nBosons; boson++)
  {
    double m2, m2prime;
    // Calculate the boson mass
    m2      = 0;
    m2prime = 0;
    BosonMasses.push_back(m2);
    // Set variables
    const double bD1      = (*Kfac)(K_type::D1, P_type::boson, sqrt(m2));
    const double bD2      = (*Kfac)(K_type::D2, P_type::boson, sqrt(m2));
    const double bR       = -Ki->vw;
    const double bRbar    = (*Kfac)(K_type::Rbar, P_type::boson, sqrt(m2));
    const double bQ1      = (*Kfac)(K_type::Q1, P_type::boson, sqrt(m2));
    const double bQ2      = (*Kfac)(K_type::Q2, P_type::boson, sqrt(m2));
    const double bm2prime = m2prime;
    // Calculate A inverse for bosons
    Ainverse[nF2 + 2 * boson][nF2 + 2 * boson]     = bR / (bD2 - bD1 * bR);
    Ainverse[nF2 + 2 * boson][nF2 + 2 * boson + 1] = -1 / (bD2 - bD1 * bR);
    Ainverse[nF2 + 2 * boson + 1][nF2 + 2 * boson] = bD2 / (bD2 - bD1 * bR);
    Ainverse[nF2 + 2 * boson + 1][nF2 + 2 * boson + 1] =
        -bD1 / (bD2 - bD1 * bR);
    // Calculate Gamma = deltaC - m2'B
    m2B[nF2 + 2 * boson][nF2 + 2 * boson] = Ki->gamw * Ki->vw * bQ1 * bm2prime;
    m2B[nF2 + 2 * boson + 1][nF2 + 2 * boson] =
        Ki->gamw * Ki->vw * bQ2 * bm2prime;
    m2B[nF2 + 2 * boson + 1][nF2 + 2 * boson + 1] = bRbar * bm2prime;
  }

  // Gamma = deltaC - m2' B
  const std::vector<std::vector<double>> CollisiontMatrix =
      CalculateCollisionMatrix(mW, FermionMasses, BosonMasses);

  // Calculate M = A^-1 * Gamma ( = deltaC - m2'B)
  for (int i = 0; i < 2 * (nBosons + nFermions); i++)
    for (int j = 0; j < 2 * (nBosons + nFermions); j++)
      for (int l = 0; l < 2 * (nBosons + nFermions); l++)
        Mtilde[i][j] += Ainverse[i][l] * (CollisiontMatrix[l][j] - m2B[l][j]);

  // Calculate Stilde = A^-1 * S
  for (int i = 0; i < 2 * (nBosons + nFermions); i++)
    for (int j = 0; j < 2 * (nFermions); j++)
      Stilde[i] += Ainverse[i][j] * S[j];
}

void TransportEquations::SolveTransportEquation()
{
  std::vector<double> zList;
  std::vector<std::vector<double>> STildeList;
  std::vector<std::vector<std::vector<double>>> MTildeList;

  std::vector<std::vector<double>> Mtilde, MtildeM1, MtildeP1;
  std::vector<double> Stilde, StildeM1, StildeP1;

  Equations(-1, MtildeM1, StildeM1);
  Equations(1, MtildeP1, StildeP1);

  for (double z = -10; z < 10; z += 1. / 1000.)
  {
    // Compute Mtilde and Stilde
    if (z < -1)
    {
      Mtilde = MtildeM1;
      Stilde = StildeM1;
    }
    else if (z > 1)
    {
      Mtilde = MtildeP1;
      Stilde = StildeP1;
    }
    else
      Equations(z, Mtilde, Stilde);

    // Save the Mtilde and Stilde
    zList.push_back(z);
    MTildeList.push_back(Mtilde);
    STildeList.push_back(Stilde);
  }
  // Construct Difeq object (S_j,n matrix)
  Difeq difeq(zList, nFermions, nBosons, MTildeList, STildeList);

  Doub itmax = 30;
  Doub conv  = 1e-10;
  Doub slowc = 1e-3;
  VecDoub scalv(zList.size(), 1);
  VecInt indexv(2 * (nFermions + nBosons));
  indexv[0] = 0;
  indexv[1] = 4;
  indexv[2] = 1;
  indexv[3] = 5;
  indexv[4] = 2;
  indexv[5] = 6;
  indexv[6] = 3;
  indexv[7] = 7;
  Int NB    = 4;
  MatDoub y(2 * (nFermions + nBosons), zList.size(), (double)0.);
  Solvde solvde(itmax, conv, slowc, scalv, indexv, NB, y, difeq);

  // Store the solution
  SolutionZ = zList;
  Solution  = y;
  // f
  PrintTransportEquation(80, "tL", "mu");
  PrintTransportEquation(80, "tR", "mu");
  PrintTransportEquation(80, "bL", "mu");
  PrintTransportEquation(80, "h", "mu");
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

  for (size_t i = 0; i < SolutionZ.value().size(); i++)
    if (abs(SolutionZ.value()[i]) < Lw.value() * multiplier)
    {
      z.push_back(SolutionZ.value()[i]);
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