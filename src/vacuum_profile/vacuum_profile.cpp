/**
 * @file
 */

#include <BSMPT/vacuum_profile/vacuum_profile.h>

namespace BSMPT
{
namespace VacuumProfileNS
{

VecInt VacuumProfile::Calcindexv()
{
  VecInt indexv(2 * dim);
  for (size_t i = 0; i < 2 * dim; i++)
  {
    // 0 1 2 3 4 5...
    indexv[i] = i;
    // switch index of field and deriv
    // 0 1 2 3 4 5-> 3 4 5 1 2 3
    indexv[i] += dim * ((i < dim) * 2 - 1) *
                 /* reordeing only necessary for dirichlet */
                 (mode == ProfileSolverMode::Field);
  }
  return indexv;
}

void VacuumProfile::LoadPath(const std::vector<double> &z_In,
                             const std::vector<std::vector<double>> &path_In)
{
  if (z_In.size() != path_In.size())
    throw("z and path have different lengths.");
  scalv = VecDoub(2 * dim, -1);
  z     = std::vector<double>(z_In.size(), 0.);
  y     = MatDoub(dim * 2, z.size(), 0.);
  for (size_t k = 0; k < path_In.size(); k++)
  {
    // load z position
    z[k] = z_In[k];
    for (size_t i = 0; i < dim; i++)
    {
      // load field
      y[dim + i][k] = path_In[k][i];
      // load field z-derivative
      if (k != 0 and k != path_In.size() - 1) // ignore first and last point
        y[i][k] = (path_In[k + 1][i] - path_In[k - 1][i]) /
                  (z_In[k + 1] - z_In[k - 1]);
      // calculate scale of quantities
      if (abs(y[i][k]) > scalv[i]) scalv[i] = abs(y[i][k]);
      if (abs(y[dim + i][k]) > scalv[dim + i])
        scalv[dim + i] = abs(y[dim + i][k]);
    }
  }
}

void VacuumProfile::CalculateProfile()
{
  size_t itmax = 10;
  double conv  = 1e-8;
  double slowc = 1;
  mode         = ProfileSolverMode::Field;

  VecInt indexv = Calcindexv();
  std::stringstream ss;

  ss << "---------------------------------------\nStarting vacuum profile "
        "calculation\n---------------------------------------\n";

  ss << "True vacuum = " << TrueVacuum << "\n";
  ss << "False vacuum = " << FalseVacuum << "\n";

  Difeq_VacuumProfile difeq_vacuumprofile(
      mode, dim, z, TrueVacuum, FalseVacuum, V, dV, Hessian);
  RelaxOde solvde(
      itmax, conv, slowc, scalv, indexv, dim, y, difeq_vacuumprofile);

  ss << "\nμ = \t" << difeq_vacuumprofile.mu << "\n";

  Logger::Write(LoggingLevel::VacuumProfile, ss.str());
}

double VacuumProfile::CalculateWidth(
    const std::vector<double> &TrueVacuum,
    const std::vector<double> &FalseVacuum,
    const std::function<double(std::vector<double>)> &V)
{
  size_t NumberPointsBarrier = 100;

  Logger::Write(LoggingLevel::VacuumProfile,
                "Estimating Lw using a rough approximation");

  const double vc     = BSMPT::L2NormVector(TrueVacuum - FalseVacuum);
  const double Vtrue  = V(TrueVacuum);
  const double Vfalse = V(FalseVacuum);
  double Vb           = -1.;
  for (size_t k = 0; k < NumberPointsBarrier; k++)
  {
    const double Vk = V(TrueVacuum + (double)k / (NumberPointsBarrier - 1.) *
                                         (FalseVacuum - TrueVacuum));
    Vb              = std::max(Vb, std::max(abs(Vtrue - Vk), abs(Vfalse - Vk)));
  }
  const double Lw = vc / sqrt(8 * Vb);

  return Lw;
}
VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t &dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In,
    // z knot list
    std::vector<double> &z_In,
    // field knot list
    std::vector<std::vector<double>> &path_In)
    : dim(dim_In)
    , TrueVacuum(TrueVacuum_In)
    , FalseVacuum(FalseVacuum_In)
    , V(V_In)
    , dV(dV_In)
    , Hessian(Hessian_In)
{
  LoadPath(z_In, path_In);
}

VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t &dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In,
    // Bubble width
    const double &Lw)
    : dim(dim_In)
    , TrueVacuum(TrueVacuum_In)
    , FalseVacuum(FalseVacuum_In)
    , V(V_In)
    , dV(dV_In)
    , Hessian(Hessian_In)
{
  // temporary variables
  std::vector<std::vector<double>> path;
  std::vector<double> vev_k(dim), zpath;
  // Use the kink to calculate a first approximation
  Logger::Write(LoggingLevel::VacuumProfile, "Lw = " + std::to_string(Lw));
  for (size_t k = 0; k < NumberOfSteps; k++)
  {
    // uniformly distributed between -10 Lw and 10 Lw
    const double zk = (-10. + (20.) * k / (NumberOfSteps - 1.)) * Lw;
    for (size_t i = 0; i < dim; i++)
      vev_k.at(i) =
          TrueVacuum.at(i) +
          (1. + tanh(zk / Lw)) / 2. * (FalseVacuum.at(i) - TrueVacuum.at(i));
    zpath.push_back(zk);
    path.push_back(vev_k);
  }
  LoadPath(zpath, path);
}

VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t &dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In)
    : VacuumProfile(dim_In,
                    TrueVacuum_In,
                    FalseVacuum_In,
                    V_In,
                    dV_In,
                    Hessian_In,
                    CalculateWidth(TrueVacuum_In, FalseVacuum_In, V_In)) {};

} // namespace VacuumProfileNS
} // namespace BSMPT