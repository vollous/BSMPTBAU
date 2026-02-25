#include <BSMPT/baryo_fhck/TransportEquations.h>

namespace BSMPT
{

namespace Baryo
{
namespace FHCK
{
TransportEquations::TransportEquations(
    const std::shared_ptr<TransportModel> &model_in,
    const double &Tstar_in,
    const std::vector<size_t> &moments_in)
    : transportmodel(model_in)
    , Tstar(Tstar_in)
    , moments(moments_in)
    , BAUeta(moments.size())
{
  stringstream ss;
  ss << "Moments to calculate = " << moments;
  Logger::Write(LoggingLevel::FHCK, ss.str());
  Initialize();
}

void TransportEquations::Initialize()
{
  transportmodel->Initialize();

  if (transportmodel->status != TransportModelStatus::Success)
  {
    Logger::Write(LoggingLevel::FHCK,
                  "Field profile calculation failed. Aborting...\n");
    return;
  }

  MakeDistribution(1, NumberOfSteps);

  Logger::Write(LoggingLevel::FHCK,
                "Limits in z \t" + std::to_string(zList.front()) + " -> " +
                    std::to_string(zList.back()));

  Logger::Write(LoggingLevel::FHCK,
                "Limits in u \t" + std::to_string(uList.front()) + " -> " +
                    std::to_string(uList.back()) + "\n");

  Logger::Write(LoggingLevel::FHCK, "Calculating fermion masses.\n");
  transportmodel->GenerateFermionMass(zList);

  gamwall = 1. / std::sqrt(1. - transportmodel->vwall * transportmodel->vwall);

  Logger::Write(LoggingLevel::FHCK, "Building Kernels Interpolations...");

  BuildKernelInterpolation();

  Logger::Write(LoggingLevel::FHCK, "\033[92mSuccess.\033[0m");
}

void TransportEquations::InitializeMoment(const size_t &moment_in)
{
  moment = moment_in;

  nParticles = nFermions + nBosons;
  nEqs       = moment * nParticles;

  Solution = MatDoub(nEqs, uList.size(), 0.);

  MtildeM.resize(nEqs, nEqs);
  MtildeP.resize(nEqs, nEqs);
  StildeM.resize(nEqs);
  StildeP.resize(nEqs);
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
        res_vals.push_back(spl(transportmodel->vwall));
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

  // Pad first element with empty spline
  Qlf.push_back(tk::spline());
  Q8ol.push_back(tk::spline());
  Q9ol.push_back(tk::spline());
  Qlb.push_back(tk::spline());

  const size_t max_moment = *std::max_element(moments.begin(), moments.end());

  for (size_t l = 0; l <= max_moment; l++)
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

MatDoub TransportEquations::CalculateCollisionMatrix(const double &mW,
                                                     VecDoub &FermionMasses,
                                                     VecDoub &BosonMasses)
{
  constexpr double KAPPABE = 0.685;
  constexpr double KAPPAFD = 0.912;
  if (FermionMasses.size() != 3)
    throw std::runtime_error(
        "Transport equation are only supported with 3 Weyl fermions.");

  if (BosonMasses.size() != 1)
    throw std::runtime_error(
        "Transport equation are only supported with 1 boson.");

  std::vector<double> mass = {sqrt(FermionMasses[0]),
                              sqrt(FermionMasses[1]),
                              sqrt(FermionMasses[2]),
                              sqrt(BosonMasses[0])};

  std::vector<ParticleType> ptype = {ParticleType::LeftFermion,
                                     ParticleType::RightFermion,
                                     ParticleType::LeftFermion,
                                     ParticleType::Boson};

  const double D0T = Dlf[0](mass[0]);
  const double D0B = Dlf[0](mass[2]);
  const double D0H = Dlb[0](mass[3]);
  const double D2T = Dlf[2](mass[0]);
  const double D2B = Dlf[2](mass[2]);
  const double D2H = Dlb[2](mass[3]);

  const std::vector<double> gammatot = {D2T / (D0T * 6.) * Tstar,
                                        D2T / (D0T * 6.) * Tstar,
                                        D2B / (D0B * 6.) * Tstar,
                                        D2H / (D0H * 20.) * Tstar};

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
  double m, Kp, ul, kappa;
  for (size_t i = 0; i < prtcl.size(); i++)
  {
    m     = mass[i];
    kappa = (ptype[i] == ParticleType::Boson ? KAPPABE : KAPPAFD);
    for (size_t l = 0; l < moment; l++)
    {
      if (l == 0)
      {
        ul                      = 0.;
        trunc[rates.size() - 1] = 1.;
        Kp = (ptype[i] == ParticleType::Boson ? Klb[l](m) : Klf[l](m));
      }
      else if (l == 1)
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 1.;
        Kp                      = -transportmodel->vwall *
             (ptype[i] == ParticleType::Boson ? Klb[0](m) : Klf[0](m));
      }
      else
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 0.;
        Kp = (ptype[i] == ParticleType::Boson ? Klb[l](m) : Klf[l](m));
      }

      for (size_t j = 0; j < prtcl.size(); j++)
      {
        for (size_t k = 0; k < rates.size(); k++)
        {
          Gamma[moment * i + l][moment * j] +=
              prtcl[i][j][k] * rates[k] * trunc[k];
        }
        Gamma[moment * i + l][moment * j] *= kappa * Kp;
      }
      Gamma[moment * i + l][moment * i + l] -= kappa * ul * gammatot[i];
    }
  }

  return Gamma;
}

std::vector<double> TransportEquations::calc_Ri(const size_t &particle,
                                                const size_t &k)
{
  std::vector<double> Ri(moment, 0);

  switch (transportmodel->truncationscheme)
  {
  case TruncationScheme::Zero:
    // 0, 0, ..., 0, -1
    Ri.at(moment - 2) = 0;
    break;
  case TruncationScheme::MinusVw:
    // 0, 0, ..., -vwall, -1
    Ri.at(moment - 2) = -transportmodel->vwall;
    break;
  case TruncationScheme::One:
    // 0, 0, ..., 1, -1
    Ri.at(moment - 2) = 1;
    break;
  case TruncationScheme::Variance:
  {
    // R_1, ..., R_n-1, -1
    double n = (double)moment; // This conversion is NECESSARY!
    Ri.at(0) = -n * n * pow(-Solution[particle * moment + 1][k], moment - 1);
    for (size_t i = 2; i < moment; i++)
    {
      double nck = nChoosek(moment, i);
      Ri.at(0) += nck * (n - i) * Solution[particle * moment + i][k] *
                  pow(-1 * Solution[particle * moment + 1][k], moment - i - 1);
      Ri.at(i - 1) =
          -nck * pow(-Solution[particle * moment + 1][k], moment - i);
    }
    break;
  }

  default: throw("Error on selecting the truncation scheme");
  }

  Ri.at(moment - 1) = -1;

  return Ri;
}

MatDoub TransportEquations::calc_Ainv(const double &m,
                                      const ParticleType &type,
                                      const size_t &particle,
                                      const size_t &k)
{
  MatDoub res(moment, moment, 0.);
  std::vector<double> Di(moment, 1), Ri = calc_Ri(particle, k);

  for (size_t i = 0; i < moment - 1; i++)
  {
    // 1, D1, D2, ..., Dn-1
    Di.at(i + 1) =
        (type == ParticleType::Boson ? Dlb[i + 1](m) : Dlf[i + 1](m));
    // off-diagonal "diagonal"
    res[i + 1][i] = 1;
  }

  // (-1)^momentInv(A)
  double Dn = (type == ParticleType::Boson ? Dlb[moment](m) : Dlf[moment](m));
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
  const double fRbar = (type == ParticleType::Boson ? Rbarb(m) : Rbarf(m));

  for (size_t l = 1; l <= moment; l++)
  {
    const double fQ   = (type == ParticleType::Boson ? Qlb[l](m) : Qlf[l](m));
    res[l - 1][l - 1] = fRbar * (double)(l - 1) * dm2;
    res[l - 1][0]     = gamwall * transportmodel->vwall * fQ * dm2;
  }
  return res;
}

VecDoub TransportEquations::calc_source(const double &m,
                                        const double &dm2,
                                        const double &dth,
                                        const double &d2th,
                                        const int &h)
{
  VecDoub res(moment);
  for (size_t i = 0; i < moment; i++)
    res[i] = -transportmodel->vwall * gamwall * h *
             ((dm2 * dth + m * m * d2th) * Q8ol[i + 1](m) -
              dm2 * m * m * dth * Q9ol[i + 1](m)); // Si
  return res;
}

void TransportEquations::InsertBlockDiagonal(MatDoub &full,
                                             MatDoub &sub,
                                             const size_t position)
{
  assert((full.cols() % sub.cols()) == 0);
  size_t start = position * sub.cols();
  for (size_t i = 0; i < sub.rows(); i++)
    for (size_t j = 0; j < sub.cols(); j++)
      full[start + i][start + j] = sub[i][j];
}

double TransportEquations::dudz(const double &u)
{
  return pow(1 - u * u, 3. / 2.) / (transportmodel->Lw * LwMultiplier);
}

double TransportEquations::zTOu(const double &z)
{
  return (z / (transportmodel->Lw * LwMultiplier)) /
         sqrt(1 + pow(z / (transportmodel->Lw * LwMultiplier), 2));
}

double TransportEquations::uTOz(const double &u)
{
  return LwMultiplier * transportmodel->Lw * u / (sqrt(1 - u * u));
}

void TransportEquations::Equations(const double &z,
                                   MatDoub &Mtilde,
                                   VecDoub &Stilde,
                                   const size_t &k)
{
  MatDoub Ainverse(nEqs, nEqs, 0.); // Store A^-1

  VecDoub S(nEqs, 0.); // Source vector

  MatDoub m2B(nEqs, nEqs, 0.); // m2' B;

  Stilde.zero(); // A^-1 * Source vector
  Mtilde.zero(); // Store A^-1 * M

  // Mass vector
  const double mW = transportmodel->GetWMass(z, Tstar);
  VecDoub FermionMasses(nFermions);
  VecDoub BosonMasses(nBosons);

  for (size_t particle = 0; particle < transportmodel->involvedparticles.size();
       particle++)
  {
    ParticleType ptype = transportmodel->involvedparticles[particle];
    double m2, m2prime, thetaprime, theta2prime;
    if (ptype == ParticleType::Boson)
    {
      m2             = 0;
      m2prime        = 0;
      BosonMasses[0] = m2;
    }
    else
    {
      // Calculate fermionic masses
      transportmodel->GetFermionMass(
          z, particle, m2, m2prime, thetaprime, theta2prime);
      FermionMasses[particle] = m2;

      // Source terms
      int h = (ptype == ParticleType::LeftFermion ? -1 : 1);
      VecDoub tempS(calc_source(sqrt(m2), m2prime, thetaprime, theta2prime, h));
      for (size_t i = 0; i < moment; i++)
        S[moment * particle + i] = tempS[i];
    }

    // Calculate A inverse for fermion
    MatDoub tempA(calc_Ainv(sqrt(m2), ptype, particle, k));
    InsertBlockDiagonal(Ainverse, tempA, particle);

    // Calculate m2'B for fermion
    MatDoub tempB(calc_m2B(sqrt(m2), m2prime, ptype));
    InsertBlockDiagonal(m2B, tempB, particle);
  }

  // Gamma = deltaC - m2' B
  const MatDoub CollisiontMatrix =
      CalculateCollisionMatrix(mW, FermionMasses, BosonMasses);

  // Calculate M = A^-1 * Gamma ( = deltaC - m2'B)
  for (size_t i = 0; i < nEqs; i++)
    for (size_t j = 0; j < nEqs; j++)
      for (size_t l = 0; l < nEqs; l++)
        Mtilde[i][j] += Ainverse[i][l] * (CollisiontMatrix[l][j] - m2B[l][j]);

  // Calculate Stilde = A^-1 * S
  for (size_t i = 0; i < nEqs; i++)
    for (size_t j = 0; j < moment * (nFermions); j++)
      Stilde[i] += Ainverse[i][j] * S[j];
}

void TransportEquations::MakeDistribution(const double xmax,
                                          const size_t npoints)
{
  uList.resize(npoints);
  for (size_t i = 0; i < npoints; i++)
  {
    double temp =
        pow(((i + 1 - (npoints + 1) / 2.)) / ((npoints + 1) / 2.), 3) * M_PI /
        4.;
    uList.at(i) = xmax * tan(temp);
  }

  zList.clear();
  for (const auto &u : uList)
  {
    zList.push_back(uTOz(u));
  }
}

void TransportEquations::CheckBoundary()
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
    ss << " \033[31m\nToo many non-decaying modes (" << NumberOfNonDecayingModes
       << ") for the boundary conditions (" << nEqs
       << "). Impossible to satisfy "
          "the "
          "boundary "
          "conditions of μ = u = 0.\033[0m\n";
    Status = FHCKStatus::UnphysicalBoundary;
  }
  else
    ss << " \033[92mSuccess.\033[0m";

  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void TransportEquations::smatrix(const int k,
                                 const int k1,
                                 const int k2,
                                 const int jsf,
                                 VecInt &indexv,
                                 MatDoub &s,
                                 MatDoub &y)
{
  double temp;
  MatDoub Mtilde(nEqs, nEqs);
  VecDoub Stilde(nEqs);
  s.zero(); // Set matrix s = 0
  if (k == k1)
  {
    // Boundary conditions mu = 0 and first u on first boundary
    for (size_t i = 0; i < nEqs / 2; i++)
    {
      // Sn at the first boundary
      s[nEqs / 2 + i][nEqs + i] = 1.0;
      // B0
      s[nEqs / 2 + i][jsf] = y[indexv[i]][0];
    }
  }
  else if (k > k2 - 1)
  {
    // Boundary conditions mu = 0 and first u on second boundary
    for (size_t i = 0; i < nEqs / 2; i++)
    {
      // Sn at the last boundary
      s[i][nEqs + i] = 1.0;
      // C0
      s[i][jsf] = y[indexv[i]][uList.size() - 1];
    }
  }
  else
  {
    double um = (uList[k] + uList[k - 1]) / 2.;
    double du = (uList[k] - uList[k - 1]) / dudz(um);
    if (um < -0.7) // TODO: change. Too specific
    {
      Mtilde = MtildeM;
      Stilde = StildeM;
    }
    else if (um > 0.7) // TODO: change. Too specific
    {
      Mtilde = MtildeP;
      Stilde = StildeP;
    }
    else
      Equations(uTOz(um), Mtilde, Stilde, k);
    for (size_t j = 0; j < nEqs; j++)
    {
      for (size_t n = 0; n < nEqs; n++)
      {
        // s matrix for the middle point
        // S_{j,n}
        // M[k] and S[k] are evaluated between k and k-1
        s[j][indexv[n]] = -Delta(j, n) - 0.5 * du * (Mtilde[j][n]);
        // S_{j,N + n}
        s[j][nEqs + indexv[n]] = Delta(j, n) - 0.5 * du * (Mtilde[j][n]);
      }
      //  Equations for E(k,k-1)
      temp = Stilde[j];
      for (size_t i = 0; i < nEqs; i++)
        temp += Mtilde[j][i] * (y[i][k] + y[i][k - 1]) / 2.;
      s[j][jsf] = y[j][k] - y[j][k - 1] - du * temp;
    }
  }
}

void TransportEquations::SolveTransportEquation()
{
  if (transportmodel->status != TransportModelStatus::Success)
  {

    Logger::Write(LoggingLevel::FHCK,
                  "SolveTransportEquation() failed. No VEV profile available.");
    std::fill(BAUeta.begin(), BAUeta.end(), NAN);
    return;
  }

  Logger::Write(LoggingLevel::FHCK,
                "Lw = " + std::to_string(transportmodel->Lw));

  for (size_t ell = 0; ell < moments.size(); ell++)
  {
    InitializeMoment(moments.at(ell));
    bau = 0;

    Equations(zList.front(), MtildeM, StildeM, 0);
    Equations(zList.back(), MtildeP, StildeP, NumberOfSteps - 1);

    CheckBoundary();
    if (Status != FHCKStatus::NotSet)
    {
      Logger::Write(
          LoggingLevel::FHCK,
          "\033[31mTransport equation unable to be computed. Status code : " +
              FHCKStatusToString.at(FHCKStatus::NotSet) + "\033[0m\n");
      return;
    }

    size_t itmax = 10;
    double conv  = 1e-10;
    double slowc = 1.;
    VecDoub scalv(nEqs, 1);
    VecInt indexv(nEqs);

    // fix μ and first u
    for (size_t l = 0; l < moment /* moment= 2 + 4k */; l++)
      for (size_t particle = 0; particle < nParticles; particle++)
      {
        indexv[l + particle * moment] = particle + l * (nParticles);
        Logger::Write(LoggingLevel::FHCK,
                      "indexv[" + std::to_string(particle + l * (nParticles)) +
                          "] = " + std::to_string(l + particle * moment));
      }

    int NB = nEqs / 2;

    RelaxOde solvde(itmax, conv, slowc, scalv, indexv, NB, Solution, *this);

    std::string str = "u.tsv";
    std::ofstream res(str);

    for (size_t k = 0; k < uList.size(); k++)
    {
      res << uList[k];
      for (size_t i = 0; i < nEqs; i++)
        res << "\t" << Solution[i][k];
      res << "\n";
    }

    res.close();

    PrintTransportEquation(120, "tL", "mu");
    PrintTransportEquation(120, "tR", "mu");
    PrintTransportEquation(120, "bL", "mu");
    PrintTransportEquation(120, "h", "mu");

    CalculateBAU();

    BAUeta.at(ell) = bau;
  }
}

void TransportEquations::CalculateBAU()
{
  double mt2, mb2, m2prime, thetaprime, theta2prime; // temporary vars
  // Weak spharelon rate
  const double Gsph = 1.e-6 * Tstar;
  double r;                // temporary variable to store the result
  std::vector<double> u;   // list of u positions
  std::vector<double> muB; // muB integrand at position u
  for (size_t i = 0; i < uList.size(); i++)
  {

    const double ui = uList[i]; // u at position i
    const double zi = uTOz(ui); // z(u)
    // Calculate the vev at z
    transportmodel->GetFermionMass(
        zi, 0, mt2, m2prime, thetaprime, theta2prime);
    transportmodel->GetFermionMass(
        zi, 2, mb2, m2prime, thetaprime, theta2prime);
    // Results
    r = 0;
    r += (1 + 4 * Dlf[0](sqrt(mt2))) / 2. * Solution[0][i];          // tL
    r += (1 + 4 * Dlf[0](sqrt(mb2))) / 2. * Solution[moment * 2][i]; // bL
    r += 2. * Dlf[0](sqrt(mt2)) * Solution[moment][i];               // tR
    r *= min(1.,
             1.7 * Tstar / Gsph *
                 exp(-37 * transportmodel->EWSBVEV(zi) / Tstar)); // f_sph(z)
    r *= exp(-45 * Gsph * std::abs(zi) /
             (4. * transportmodel->vwall *
              gamwall));    // exp(-45 G_sph |z| / 4 vw gammaw)
    r *= 1 / abs(dudz(ui)); // z -> u jacobian
    // Save in list to pass to integrator
    u.push_back(ui);
    muB.push_back(r); // integrand
  }
  // Step 1: Set up GSL interpolation
  size_t n = u.size();
  const gsl_interp_type *interp_type =
      gsl_interp_cspline; // Cubic spline interpolation
  gsl_interp *interp = gsl_interp_alloc(interp_type, n);
  gsl_interp_init(interp, u.data(), muB.data(), n);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  // Step 2: Integrate the interpolated function
  double a = u.front(); // Lower bound of integration
  double b = u.back();  // Upper bound of integration
  double result, error;
  // Struct to hold parameters
  struct Params
  {
    gsl_interp *interp;
    gsl_interp_accel *acc;
    const std::vector<double> *u;
    const std::vector<double> *muB;
  } params                             = {interp, acc, &u, &muB};
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = [](double x, void *p) -> double
  {
    auto *data = static_cast<Params *>(p);
    return gsl_interp_eval(
        data->interp, data->u->data(), data->muB->data(), x, data->acc);
  };
  F.params = &params;
  gsl_integration_qags(
      &F, a, b, 1e-30, 1e-30, 1000, workspace, &result, &error);
  // Cleanup
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  // Results
  const double prefactor = 405. * Gsph /
                           (4 * M_PI * M_PI * transportmodel->vwall * gamwall *
                            106.75 * Tstar); // TODO Fix gstas

  result *= prefactor;
  error *= prefactor;
  // Save the result
  bau = result;
  // Step 3: Output the result
  stringstream ss;
  ss << "eta = " << result << " with error " << error << std::endl;

  double unc = abs(error / result);

  if (unc > 0.01)
  {
    ss << "Calculation failed!\t" << unc << "\n";
    result = NAN;
  }

  ss << "eta/eta_obs = " << result / (8.7e-11) << std::endl;
  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void TransportEquations::PrintTransportEquation(const int &size,
                                                const std::string &Particle,
                                                const std::string &MuOrU)
{
  AsciiPlotter Plot(Particle + " " + MuOrU, size, ceil(size / 3.));
  std::optional<int> ind;
  std::vector<double> z, y;

  if (Particle == "tL") ind = 0 * moment;
  if (Particle == "tR") ind = 1 * moment;
  if (Particle == "bL") ind = 2 * moment;
  if (Particle == "h") ind = 3 * moment;

  if (not ind.has_value()) throw("Invalid particle to plot the solution.");

  if (MuOrU == "u") ind = ind.value() + 1;

  size_t i_min_left = 0, i_min_right = uList.size() - 1;
  double max = -1.;
  for (size_t i = 0; i < uList.size(); i++)
    max = std::max(max, abs(Solution[ind.value()][i]));

  for (size_t i = 0; i < uList.size(); i++)
  {
    if (abs(Solution[ind.value()][i]) > max / 100.)
    {
      i_min_left = i;
      break;
    }
  }

  for (size_t i = uList.size() - 1; i >= 0; i--)
  {
    if (abs(Solution[ind.value()][i]) > max / 100.)
    {
      i_min_right = i;
      break;
    }
  }

  for (size_t i = i_min_left; i <= i_min_right; i++)
  {
    z.push_back(uList[i]);
    y.push_back(Solution[ind.value()][i]);
  }
  Plot.addPlot(z, y, "", '*');
  std::stringstream ss;
  Plot.show(ss);
  Logger::Write(LoggingLevel::FHCK, ss.str());
}
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT
