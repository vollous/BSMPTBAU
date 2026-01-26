#include <BSMPT/baryo_fhck/TransportEquations.h>

namespace BSMPT
{

namespace Baryo
{
namespace FHCK
{
TransportEquations::TransportEquations(
    const std::shared_ptr<TransportModel> &model_in,
    const double &Tstar_in)
{
  Tmodel = model_in;
  Tstar  = Tstar_in;
  Initialize();
}

void TransportEquations::Initialize()
{
  uList = MakeDistribution(1, NumberOfSteps);

  zList.clear();
  for (const auto &u : uList)
  {
    zList.push_back(uTOz(u));
  }

  Logger::Write(LoggingLevel::FHCK,
                "Limits in z \t" + std::to_string(zList.front()) + " -> " +
                    std::to_string(zList.back()) + "\n");

  Logger::Write(LoggingLevel::FHCK,
                "Limits in u \t" + std::to_string(uList.front()) + " -> " +
                    std::to_string(uList.back()) + "\n");

  Logger::Write(LoggingLevel::FHCK, "Calculating fermion masses.\n");

  gamwall = 1. / std::sqrt(1. - Tmodel->vwall * Tmodel->vwall);

  nEqs = moment * (nFermions + nBosons);

  Logger::Write(LoggingLevel::FHCK, "Building Kernels Interpolations\n");

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
        res_vals.push_back(spl(Tmodel->vwall));
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

  std::vector<ParticleType> ptype = {ParticleType::Fermion,
                                     ParticleType::Fermion,
                                     ParticleType::Fermion,
                                     ParticleType::Boson};

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
        Kp =
            (ptype[i] == ParticleType::Fermion ? Klf[l - 1](m) : Klb[l - 1](m));
      }
      else if (l == 2)
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 1.;
        Kp                      = -Tmodel->vwall *
             (ptype[i] == ParticleType::Fermion ? Klf[0](m) : Klb[0](m));
      }
      else
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 0.;
        Kp =
            (ptype[i] == ParticleType::Fermion ? Klf[l - 1](m) : Klb[l - 1](m));
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
    Di.at(i + 1) =
        (type == ParticleType::Fermion ? Dlf[i + 1](m) : Dlb[i + 1](m));
    // off-diagonal "diagonal"
    res[i + 1][i] = 1;
  }

  // 0, 0, ..., -vwall, -1
  Ri.at(moment - 2) = -Tmodel->vwall;
  Ri.at(moment - 1) = -1;

  // (-1)^moment Inv(A)
  double Dn = (type == ParticleType::Fermion ? Dlf[moment](m) : Dlb[moment](m));
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
  const double fRbar = (type == ParticleType::Fermion ? Rbarf(m) : Rbarb(m));

  for (size_t l = 1; l <= moment; l++)
  {
    const double fQ   = (type == ParticleType::Fermion ? Qlf[l](m) : Qlb[l](m));
    res[l - 1][l - 1] = fRbar * (double)(l - 1) * dm2;
    res[l - 1][0]     = gamwall * Tmodel->vwall * fQ * dm2;
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
    res[i] = -Tmodel->vwall * gamwall * h *
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
  return pow(1 - u * u, 3. / 2.) / (Tmodel->Lw * LwMultiplier);
}

double TransportEquations::zTOu(const double &z)
{
  return (z / (Tmodel->Lw * LwMultiplier)) /
         sqrt(1 + pow(z / (Tmodel->Lw * LwMultiplier), 2));
}

double TransportEquations::uTOz(const double &u)
{
  return LwMultiplier * Tmodel->Lw * u / (sqrt(1 - u * u));
}

void TransportEquations::Equations(const double &z,
                                   MatDoub &Mtilde,
                                   VecDoub &Stilde)
{
  MatDoub Ainverse(nEqs, nEqs, 0.); // Store A^-1

  VecDoub S(nEqs, 0.); // Source vector

  MatDoub m2B(nEqs, nEqs, 0.); // m2' B;

  Stilde.zero(); // A^-1 * Source vector
  Mtilde.zero(); // Store A^-1 * M

  // Mass vector
  const double mW = Tmodel->GetWMass(z, Tstar);
  VecDoub FermionMasses(nFermions);
  VecDoub BosonMasses(nBosons);

  for (size_t particle = 0; particle < Tmodel->Tprtcls.size(); particle++)
  {
    ParticleType ptype = (Tmodel->Tprtcls[particle] == TransportParticles::Boson
                              ? ParticleType::Boson
                              : ParticleType::Fermion);
    double m2, m2prime, thetaprime, theta2prime;
    if (ptype == ParticleType::Fermion)
    {
      // Calculate fermionic masses
      Tmodel->GetFermionMass(z, particle, m2, m2prime, thetaprime, theta2prime);
      FermionMasses[particle] = m2;

      // Source terms
      int h =
          (Tmodel->Tprtcls[particle] == TransportParticles::LeftFermion ? -1
                                                                        : 1);
      VecDoub tempS(calc_source(sqrt(m2), m2prime, thetaprime, theta2prime, h));
      for (size_t i = 0; i < moment; i++)
        S[moment * particle + i] = tempS[i];
    }
    else
    {
      m2             = 0;
      m2prime        = 0;
      BosonMasses[0] = m2;
    }

    // Calculate A inverse for fermion
    MatDoub tempA(calc_Ainv(sqrt(m2), ptype));
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

std::vector<double> TransportEquations::MakeDistribution(const double xmax,
                                                         const size_t npoints)
{
  std::vector<double> res(npoints);

  for (size_t i = 0; i < npoints; i++)
  {
    double temp =
        pow(((i + 1 - (npoints + 1) / 2.)) / ((npoints + 1) / 2.), 3) * M_PI /
        4.;
    res.at(i) = xmax * tan(temp);
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

void TransportEquations::SolveTransportEquation()
{
  Logger::Write(LoggingLevel::FHCK, "Lw = " + std::to_string(Tmodel->Lw));

  MatDoub STildeList(NumberOfSteps, nEqs);
  Mat3DDoub MTildeList(NumberOfSteps, nEqs, nEqs);

  MatDoub Mtilde(nEqs, nEqs), MtildeM(nEqs, nEqs), MtildeP(nEqs, nEqs);
  VecDoub Stilde(nEqs), StildeM(nEqs), StildeP(nEqs);

  Equations(zList.front(), MtildeM, StildeM);
  Equations(zList.back(), MtildeP, StildeP);

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
    double uc = (uList[i] + uList[i - 1]) / 2.;
    if (uc < -0.3) // TODO: Fix this, too specific.
    {
      Mtilde = MtildeM;
      Stilde = StildeM;
    }
    else if (uc > 0.3)
    {
      Mtilde = MtildeP;
      Stilde = StildeP;
    }
    else
      Equations(uTOz(uc), Mtilde, Stilde);

    // Save the Mtilde and Stilde
    for (size_t j = 0; j < nEqs; j++)
    {
      STildeList[i][j] = Stilde[j] / dudz(uc);
      for (size_t k = 0; k < nEqs; k++)
        MTildeList[i][j][k] = Mtilde[j][k] / dudz(uc);
    }
  }
  // Construct Difeq object (S_j,n matrix)
  Difeq_TransportEquation difeq(
      uList, nFermions, nBosons, MTildeList, STildeList, moment);

  size_t itmax = 1;
  double conv  = 1e-10;
  double slowc = 1.;
  VecDoub scalv(nEqs, 1);
  VecInt indexv(nEqs);

  // fix μ and first u
  for (size_t l = 0; l < moment /* moment = 2 + 4k */; l++)
    for (size_t particle = 0; particle < nFermions + nBosons; particle++)
    {
      indexv[l + particle * moment] = particle + l * (nFermions + nBosons);
      Logger::Write(LoggingLevel::FHCK,
                    "indexv[" +
                        std::to_string(particle + l * (nFermions + nBosons)) +
                        "] = " + std::to_string(l + particle * moment));
    }

  int NB = nEqs / 2;
  MatDoub y(nEqs, uList.size(), 0.);
  RelaxOde solvde(itmax, conv, slowc, scalv, indexv, NB, y, difeq);

  std::string str = "u.tsv";
  std::ofstream res(str);

  for (size_t k = 0; k < uList.size(); k++)
  {
    res << uList[k];
    for (size_t i = 0; i < nEqs; i++)
      res << "\t" << y[i][k];
    res << "\n";
  }

  res.close();

  // Store the solution
  Solution = y;

  PrintTransportEquation(120, "tL", "mu");
  PrintTransportEquation(120, "tR", "mu");
  PrintTransportEquation(120, "bL", "mu");
  PrintTransportEquation(120, "h", "mu");

  CalculateBAU();
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
    Tmodel->GetFermionMass(zi, 0, mt2, m2prime, thetaprime, theta2prime);
    Tmodel->GetFermionMass(zi, 2, mb2, m2prime, thetaprime, theta2prime);
    // Results
    r = 0;
    r += (1 + 4 * Dlf[0](sqrt(mt2))) / 2. * Solution.value()[0][i]; // tL
    r += (1 + 4 * Dlf[0](sqrt(mb2))) / 2. *
         Solution.value()[moment * 2][i];                      // bL
    r += 2. * Dlf[0](sqrt(mt2)) * Solution.value()[moment][i]; // tR
    r *= min(1.,
             2.4 * Tstar / Gsph *
                 exp(-40 * Tmodel->EWSBVEV(zi) / Tstar)); // f_sph(z)
    r *=
        exp(-45 * Gsph * std::abs(zi) /
            (4. * Tmodel->vwall * gamwall)); // exp(-45 G_sph |z| / 4 vw gammaw)
    r *= 1 / abs(dudz(ui));                  // z -> u jacobian
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
                           (4 * M_PI * M_PI * Tmodel->vwall * gamwall * 106.75 *
                            Tstar); // TODO Fix gstas

  result *= prefactor;
  // Save the result
  BAUEta = result;
  // Step 3: Output the result
  stringstream ss;
  ss << "eta = " << result << " with error " << error << std::endl;
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

  for (size_t i = 0; i < uList.size(); i++)
  {
    z.push_back(uList[i]);
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