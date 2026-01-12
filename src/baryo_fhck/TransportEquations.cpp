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
    const VevProfileMode &Mode_in)
{
  modelPointer = pointer_in;
  Bounce       = Bounce_in;
  Tstar        = Bounce->GetTransitionTemp();
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
  nFB2 = 2 * (nFermions + nBosons);
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
  m2      = std::abs(m * m) / (Tstar * Tstar); // m^2 of fermion
  m2prime = std::abs(2. * mprime * m) / (Tstar * Tstar);
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

  const double mTop   = sqrt(FermionMasses[0]);
  const double mBot   = sqrt(FermionMasses[2]);
  const double mHiggs = sqrt(BosonMasses[0]);

  const double D0t = (*Kfac)(K_type::D0, P_type::fermion, mTop);
  const double D0b = (*Kfac)(K_type::D0, P_type::fermion, mBot);

  const double K0t = (*Kfac)(K_type::K0, P_type::fermion, mTop); // TODO: fix
  const double K0b = (*Kfac)(K_type::K0, P_type::fermion, mBot); // TODO: fix
  const double K0h = (*Kfac)(K_type::K0, P_type::boson, mHiggs); // TODO: fix

  const double GammaTotTop = (*Kfac)(K_type::K4, P_type::fermion, mTop) /
                             ((*Kfac)(K_type::D0, P_type::fermion, mTop) * 6.) *
                             Tstar; // TODO: fix
  const double GammaTotBot = (*Kfac)(K_type::K4, P_type::fermion, mBot) /
                             ((*Kfac)(K_type::D0, P_type::fermion, mBot) * 6.) *
                             Tstar; // TODO: fix
  const double GammaTotHiggs =
      (*Kfac)(K_type::K4, P_type::boson, mHiggs) /
      ((*Kfac)(K_type::D0, P_type::boson, mHiggs) * 20.) * Tstar; // TODO: fix

  const double GammaM       = pow(mTop, 2) / 63. * Tstar;
  const double GammaY       = 4.2e-3 * Tstar;
  const double GammaW       = GammaTotHiggs;
  const double GammaTildeSS = 4.9e-4 * Tstar;
  // const double GammaTildeY  = 4.9e-4 * Tstar; Not sure what this is
  const double GammaH = mW * mW / 50. * Tstar; // TODO: fix

  MatDoub Gamma(8, 8, 0.);
  Gamma[0][0] = (GammaM + GammaW + GammaY + (1 + 9 * D0t) * GammaTildeSS) * K0t;
  Gamma[0][2] = (-GammaM - GammaY + (-1 + 9 * D0t) * GammaTildeSS) * K0t;
  Gamma[0][4] = (-GammaW + (1 + 9 * D0b) * GammaTildeSS) * K0t;
  Gamma[0][6] = GammaY * K0t;

  Gamma[1][0] = -Ki->vw * Gamma[0][0];
  Gamma[1][1] = -GammaTotTop;
  Gamma[1][2] = -Ki->vw * Gamma[0][2];
  Gamma[1][4] = -Ki->vw * Gamma[0][4];
  Gamma[1][6] = -Ki->vw * Gamma[0][6];

  Gamma[2][0] = (-GammaM - GammaY - (1 + 9 * D0t) * GammaTildeSS) * K0t;
  Gamma[2][2] = (GammaM + 2 * GammaY + (1 - 9 * D0t) * GammaTildeSS) * K0t;
  Gamma[2][4] = (-GammaY - (1 + 9 * D0b) * GammaTildeSS) * K0t;
  Gamma[2][6] = -2 * GammaY * K0t;

  Gamma[3][0] = -Ki->vw * Gamma[2][0];
  Gamma[3][2] = -Ki->vw * Gamma[2][2];
  Gamma[3][3] = -GammaTotTop;
  Gamma[3][4] = -Ki->vw * Gamma[2][4];
  Gamma[3][6] = -Ki->vw * Gamma[2][6];

  Gamma[4][0] = (-GammaW + (1 + 9 * D0t) * GammaTildeSS) * K0b;
  Gamma[4][2] = (-GammaY + (-1 + 9 * D0t) * GammaTildeSS) * K0b;
  Gamma[4][4] = (GammaW + GammaY + (1 + 9 * D0b) * GammaTildeSS) * K0b;
  Gamma[4][6] = GammaY * K0b;

  Gamma[5][0] = -Ki->vw * Gamma[4][0];
  Gamma[5][2] = -Ki->vw * Gamma[4][2];
  Gamma[5][4] = -Ki->vw * Gamma[4][4];
  Gamma[5][5] = -GammaTotBot;
  Gamma[5][6] = -Ki->vw * Gamma[4][6];

  Gamma[6][0] = GammaY * K0h;
  Gamma[6][2] = -2 * GammaY * K0h;
  Gamma[6][4] = GammaY * K0h;
  Gamma[6][6] = (GammaH + 2 * GammaY) * K0h;

  Gamma[7][0] = -Ki->vw * Gamma[6][0];
  Gamma[7][2] = -Ki->vw * Gamma[6][2];
  Gamma[7][4] = -Ki->vw * Gamma[6][4];
  Gamma[7][6] = -Ki->vw * Gamma[6][6];
  Gamma[7][7] = -GammaTotHiggs;

  return Gamma;
}

MatDoub TransportEquations::calc_Ainv(const double &m, const P_type &type)
{
  MatDoub res(2, 2);
  const double fD1 = (*Kfac)(K_type::D1, type, m);
  const double fD2 = (*Kfac)(K_type::D2, type, m);
  const double fR  = -Ki->vw;
  res[0][0]        = fR / (fD2 - fD1 * fR);
  res[0][1]        = -1 / (fD2 - fD1 * fR);
  res[1][0]        = fD2 / (fD2 - fD1 * fR);
  res[1][1]        = -fD1 / (fD2 - fD1 * fR);
  return res;
}

MatDoub TransportEquations::calc_m2B(const double &m,
                                     const double &dm2,
                                     const P_type &type)
{
  MatDoub res(2, 2);
  const double fRbar = (*Kfac)(K_type::Rbar, type, m);
  const double fQ1   = (*Kfac)(K_type::Q1, type, m);
  const double fQ2   = (*Kfac)(K_type::Q2, type, m);
  res[0][0]          = Ki->gamw * Ki->vw * fQ1 * dm2;
  res[1][0]          = Ki->gamw * Ki->vw * fQ2 * dm2;
  res[1][1]          = fRbar * dm2;
  return res;
}

VecDoub TransportEquations::calc_source(const double &m,
                                        const double &dm2,
                                        const double &dth,
                                        const double &d2th,
                                        const P_type &type)
{
  VecDoub res(2);
  const double fQ8o1 = (*Kfac)(K_type::Q8o1, type, m);
  const double fQ8o2 = (*Kfac)(K_type::Q8o2, type, m);
  const double fQ9o1 = (*Kfac)(K_type::Q9o1, type, m);
  const double fQ9o2 = (*Kfac)(K_type::Q9o2, type, m);

  res[0] =
      -Ki->vw * Ki->gamw *
      ((dm2 * dth + m * m * d2th) * fQ8o1 - dm2 * m * m * dth * fQ9o1); // S1
  res[1] =
      -Ki->vw * Ki->gamw *
      ((dm2 * dth + m * m * d2th) * fQ8o2 - dm2 * m * m * dth * fQ9o2); // S2
  return res;
}

void TransportEquations::Equations(const double &z,
                                   MatDoub &Mtilde,
                                   VecDoub &Stilde)
{
  const int nF2 = 2 * nFermions; // Clearer code

  MatDoub Ainverse(nFB2, nFB2, 0.); // Store A^-1

  VecDoub S(nFB2, 0.); // Source vector

  MatDoub m2B(nFB2, nFB2, 0.); // m2' B;

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
  for (int fermion = 0; fermion < nFermions; fermion++)
  {

    double m2, m2prime, thetaprime, theta2prime;
    // Calculate fermionic masses
    GetFermionMass(z, fermion, esQuark, m2, m2prime, thetaprime, theta2prime);
    FermionMasses[fermion] = m2;

    // Calculate A inverse for fermion
    MatDoub tempA(calc_Ainv(sqrt(m2), P_type::fermion));
    Ainverse[2 * fermion][2 * fermion]         = tempA[0][0];
    Ainverse[2 * fermion][2 * fermion + 1]     = tempA[0][1];
    Ainverse[2 * fermion + 1][2 * fermion]     = tempA[1][0];
    Ainverse[2 * fermion + 1][2 * fermion + 1] = tempA[1][1];

    // Calculate m2'B for fermion
    MatDoub tempB(calc_m2B(sqrt(m2), m2prime, P_type::fermion));
    m2B[2 * fermion][2 * fermion]         = tempB[0][0];
    m2B[2 * fermion + 1][2 * fermion]     = tempB[1][0];
    m2B[2 * fermion + 1][2 * fermion + 1] = tempB[1][1];

    // Source terms
    VecDoub tempS(calc_source(
        sqrt(m2), m2prime, thetaprime, theta2prime, P_type::fermion));
    S[2 * fermion]     = tempS[0];
    S[2 * fermion + 1] = tempS[1];
  }

  // Bosons
  for (int boson = 0; boson < nBosons; boson++)
  {
    double m2, m2prime;
    // Calculate the boson mass
    m2                 = 0;
    m2prime            = 0;
    BosonMasses[boson] = m2;

    // Calculate A inverse for bosons
    MatDoub tempA(calc_Ainv(sqrt(m2), P_type::boson));
    Ainverse[nF2 + 2 * boson][nF2 + 2 * boson]         = tempA[0][0];
    Ainverse[nF2 + 2 * boson][nF2 + 2 * boson + 1]     = tempA[0][1];
    Ainverse[nF2 + 2 * boson + 1][nF2 + 2 * boson]     = tempA[1][0];
    Ainverse[nF2 + 2 * boson + 1][nF2 + 2 * boson + 1] = tempA[1][1];

    // Calculate m2'B
    MatDoub tempB(calc_m2B(sqrt(m2), m2prime, P_type::fermion));
    m2B[nF2 + 2 * boson][nF2 + 2 * boson]         = tempB[0][0];
    m2B[nF2 + 2 * boson + 1][nF2 + 2 * boson]     = tempB[1][0];
    m2B[nF2 + 2 * boson + 1][nF2 + 2 * boson + 1] = tempB[1][1];
  }

  // Gamma = deltaC - m2' B
  const MatDoub CollisiontMatrix =
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

VecDoub TransportEquations::MakeDistribution(const double xmax,
                                             const size_t npoints)
{
  VecDoub res(npoints);

  for (size_t i = 0; i < npoints; i++)
  {
    double temp = pow(((i - npoints / 2.)) / (npoints / 2.), 3) * M_PI / 4.;
    res[i]      = xmax * tan(temp);
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

  Eigen::MatrixXcd EigenMtildeM(nFB2, nFB2);
  Eigen::MatrixXcd EigenMtildeP(nFB2, nFB2);

  for (size_t i = 0; i < nFB2; i++)
  {
    STildeLength += pow(StildeM[i], 2);
    STildeLength += pow(StildeP[i], 2);
    for (size_t j = 0; j < nFB2; j++)
    {
      EigenMtildeM(i, j) = MtildeM[i][j];
      EigenMtildeP(i, j) = MtildeP[i][j];
    }
  }

  // Checking if S vector are small enough
  STildeLength = sqrt(STildeLength);
  STildeLength /= 2. * nFB2;

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

  if (NumberOfNonDecayingModes > nFB2)
  {
    ss << " \033[31m\nToo many non-decaying modes. Impossible to satisfy the "
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
  VecDoub zList(MakeDistribution(LwMultiplier * Lw.value(), NumberOfSteps));

  Logger::Write(LoggingLevel::FHCK, "Lw = " + std::to_string(Lw.value()));

  MatDoub STildeList(NumberOfSteps, nFB2);
  Mat3DDoub MTildeList(NumberOfSteps, nFB2, nFB2);

  MatDoub Mtilde(nFB2, nFB2), MtildeM(nFB2, nFB2), MtildeP(nFB2, nFB2);
  VecDoub Stilde(nFB2), StildeM(nFB2), StildeP(nFB2);

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
    for (size_t j = 0; j < nFB2; j++)
    {
      STildeList[i][j] = Stilde[j];
      for (size_t k = 0; k < nFB2; k++)
        MTildeList[i][j][k] = Mtilde[j][k];
    }
  }
  // Construct Difeq object (S_j,n matrix)
  Difeq_TransportEquation difeq(
      zList, nFermions, nBosons, MTildeList, STildeList);

  double itmax = 30;
  double conv  = 1e-10;
  double slowc = 1e-3;
  VecDoub scalv(zList.size(), 1);
  VecInt indexv(nFB2);
  indexv[0] = 0;
  indexv[1] = 4;
  indexv[2] = 1;
  indexv[3] = 5;
  indexv[4] = 2;
  indexv[5] = 6;
  indexv[6] = 3;
  indexv[7] = 7;
  int NB    = 4;
  MatDoub y(nFB2, zList.size(), 0.);
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
  SolutionZ = zList;
  Solution  = y;

  PrintTransportEquation(120, "tL", "mu", 100);
  PrintTransportEquation(120, "tR", "mu", 100);
  PrintTransportEquation(120, "bL", "mu", 100);
  PrintTransportEquation(120, "h", "mu", 100);

  CalculateBAU();
}

void TransportEquations::CalculateBAU()
{
  double m2, m2prime, thetaprime, theta2prime; // temporary vars
  // Weak spharelon rate
  const double Gsph = 1.e-6 * Tstar;
  double r;                // temporary variable to store the result
  std::vector<double> z;   // list of z positions
  std::vector<double> muB; // muB integrand at position z
  for (size_t i = 0; i < SolutionZ.value().size(); i++)
  {

    const double zi = SolutionZ.value()[i]; // z at position i
    // Calculate the vev at z
    const std::vector<double> vev = Vev(zi, 0);
    // Calculate masses
    Eigen::MatrixXcd MIJQuarks =
        modelPointer->QuarkMassMatrix(modelPointer->MinimizeOrderVEV(vev));
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esQuark(MIJQuarks);
    // D0 of t
    GetFermionMass(zi, 0, esQuark, m2, m2prime, thetaprime, theta2prime);
    const double D0t = (*Kfac)(K_type::D0, P_type::fermion, sqrt(m2));
    // D0 of b
    GetFermionMass(zi, 2, esQuark, m2, m2prime, thetaprime, theta2prime);
    const double D0b = (*Kfac)(K_type::D0, P_type::fermion, sqrt(m2));
    // Results
    r = 0;
    r += (1 + 4 * D0t) / 2. * Solution.value()[0][i]; // tL
    r += (1 + 4 * D0b) / 2. * Solution.value()[4][i]; // bL
    r += 2. * D0b * Solution.value()[2][i];           // tR
    r *= min(
        1.,
        2.4 * Tstar / Gsph *
            exp(-40 *
                modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(vev), 0.) /
                Tstar)); // f_sph(z)
    r *= exp(-45 * Gsph * std::abs(zi) /
             (4. * Ki->vw * Ki->gamw)); // exp(-45 G_sph |z| / 4 vw gammaw)
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
      (4 * M_PI * M_PI * Ki->vw * Ki->gamw * 106.75 * Tstar); // TODO Fix gstas

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