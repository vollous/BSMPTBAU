#pragma once

#include <BSMPT/Kfactors/Kernels.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/baryo_fhck/difeq_transport_equations.h>
#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/bounce_solution/bounce_solution.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/relaxation/solvde.h>
#include <BSMPT/utility/spline/spline.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/vacuum_profile.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{
/**
 * @brief Types of vev profiles
 *
 */
enum class VevProfileMode
{
  Unset,
  Kink,
  FieldEquation
};

/**
 * @brief Status of the FHCK baryo calculation
 *
 */
enum class FHCKStatus
{
  NotSet,
  Success,
  SmallIntegrationRegion,
  UnphysicalBoundary
};

/**
 * @brief Convert FHCK baryo into strings
 *
 */
const std::unordered_map<FHCKStatus, std::string> FHCKStatusToString{
    {FHCKStatus::NotSet, "not_set"},
    {FHCKStatus::Success, "success"},
    {FHCKStatus::SmallIntegrationRegion, "small_integration_region"},
    {FHCKStatus::UnphysicalBoundary, "unphysical_boundary"}};

class TransportEquations
{
public:
  /**
   * @brief \f$ \eta = \frac{n_B}{n_\gamma}\f$
   *
   */
  std::optional<double> BAUEta;

  /**
   * @brief Status of the FHCK baryo calculation
   *
   */
  FHCKStatus Status = FHCKStatus::NotSet;

  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief pointer to the vacuum profile solver
   *
   */
  std::shared_ptr<VacuumProfileNS::VacuumProfile> vacuumprofile;

  /**
   * @brief Transition temperature
   *
   */
  double Tstar;

  /**
   * @brief Transition temperature
   *
   */
  const size_t moment = 2;

  /**
   * @brief Interpolated kernel functions for different moments
   *
   */
  std::vector<tk::spline> Klf, Klb, Dlf, Dlb, Qlf, Qlb, Q8ol, Q9ol;

  /**
   * @brief Special interpolated kernel functions
   *
   */
  tk::spline Rbarb, Rbarf, K4FHf, K4FHb;

  /**
   * @brief Bubble wall velocity
   *
   */
  double vwall;

  /**
   * @brief Rel. gamma factor of the wall
   *
   */
  double gamwall;

  /**
   * @brief Wall thickness
   *
   */
  std::optional<double> Lw;

  /**
   * @brief List of point on the z-axis.
   *
   */
  std::vector<double> zList;

  /**
   * @brief EtaInterface object. To  use BSMPTv2 functions
   *
   */
  std::shared_ptr<CalculateEtaInterface> EtaInterface;

  /**
   * @brief Phase of the top mass at the false vacuum
   *
   */
  std::optional<double> Theta_False;

  /**
   * @brief Phase of the top mass at the true vacuum
   *
   */
  std::optional<double> Theta_True;

  /**
   * @brief Store the solution from the relaxation method
   *
   */
  std::optional<MatDoub> Solution;

  /**
   * @brief Number of Weil fermions
   *
   */
  size_t nFermions = 3;

  /**
   * @brief Number of bosons
   *
   */
  size_t nBosons = 1;

  /**
   * @brief Number of equations
   *
   */
  size_t nEqs;

  /**
   * @brief Number of steps in space
   *
   */
  size_t NumberOfSteps = 2000;

  /**
   * @brief The integration goes from \f$ - LwMultiplier * Lw \f$ up to \f$
   * LwMultiplier * Lw \f$
   *
   */
  double LwMultiplier = 100.;

  /**
   * @brief Threshold for which the length of the S vector must be smaller. If
   * not then the integration region must be inscreased.
   *
   */
  double STildeThreshold = 1e-10;

  /**
   * @brief False vacuum
   *
   */
  std::vector<double> FalseVacuum;

  /**
   * @brief True vacuum
   *
   */
  std::vector<double> TrueVacuum;

  /**
   * @brief Empty vacuum
   *
   */
  std::vector<double> EmptyVacuum;

  /**
   * @brief Mode of the vev profile.
   * - Kink -> Kink solution
   * - FieldEquation -> Calculate the vev profile using the tunneling solution
   */
  VevProfileMode VevProfile = VevProfileMode::Unset;

  /**
   * @brief Coex phase object
   *
   */
  std::shared_ptr<CoexPhases> CoexPhase;

  /**
   * @brief Splines to store the variation of the quark masses
   *
   */
  std::vector<tk::spline> QuarkMassesRe, QuarkMassesIm;

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param pointer_in Model pointer
   * @param CoexPhase_in Coexphase pointer
   * @param vwall_in Wall velocity
   * @param Tstar_in Transition temperature
   * @param VevProfile_In Solver mode. Default: kink solution
   */
  TransportEquations(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const std::shared_ptr<CoexPhases> &CoexPhase_in,
      const double &vwall_in,
      const double &Tstar_in,
      const VevProfileMode &VevProfile_In = VevProfileMode::Kink);

  /**
   * @brief Create the VEV vectors and initalize **Ki()** and **Kfac()**
   * objects.
   *
   */
  void Initialize();

  tk::spline InterpolateKernel(const std::string &kernel_Name,
                               const bool is_1D);

  /**
   * @brief Load all the Kernel functions and interpolate them
   *
   *
   */
  void BuildKernelInterpolation();

  /**
   * @brief Set the number of intergrations steps and initialize again.
   *
   * @param Num number of steps
   */
  void SetNumberOfSteps(const int &num);

  /**
   * @brief Set EtaInterface obejct to use BSMPTv2 functions.
   *
   */
  void SetEtaInterface();

  /**
   * @brief Calculates the VEV \f$ \vec{v} \f$ and derivatives \f$
   * \drac{d^n\vec{v}}{dz^n} \f$ at the positions z.
   *
   * @param z
   * @param diff 0 = vev, 1 = \f$ \drac{d\vec{v}}{dz} \f$, 2 = \f$
   * \drac{d^2\vec{v}}{dz^2} \f$
   * @return std::vector<double> result
   */
  std::vector<double> Vev(const double &z, const int &diff = 0);

  /**
   * @brief Generate the quark masses splines
   *
   */
  void GenerateFermionMass();

  /**
   * @brief Get the Fermion Mass object Calculate the fermion mass and its
   * derivatives
   *
   * @param z distance to the bubble wall
   * @param fermion which fermion, 0 = most massive
   * @param m2 \f$ m^2 \f$
   * @param m2prime \f$ m'^2 \f$
   * @param thetaprime \f$ \theta' \f$
   * @param theta2prime \f$ \theta'' \f$
   */
  void GetFermionMass(const double &z,
                      const size_t &fermion,
                      double &m2,
                      double &m2prime,
                      double &thetaprime,
                      double &theta2prime);

  /**
   * @brief Calculate the W boson mass. (same code as BSMPTv2)
   *
   * @param vev VEV
   * @param T Transition temperature
   * @return double W boson mass
   */
  double GetWMass(const std::vector<double> &vev, const double &T) const;

  /**
   * @brief Calculate the collision matrix \f$ \deltaC^{1/2} \f$
   *
   * @param z distance to the bubble wall.
   * @return MatDoub Collision matrix
   */
  MatDoub CalculateCollisionMatrix(const double &mW,
                                   VecDoub &FermionMasses,
                                   VecDoub &BosonMasses);
  /**
   * @brief Calculate the 2x2 submatrix of Ainv for 1 particle
   *
   * @param m mass of the particle
   * @param type type of particle e.g. fermion/boson
   * @return MatDoub of Ainv
   */
  MatDoub calc_Ainv(const double &m, const ParticleType &type);

  /**
   * @brief Calculate the 2x2 submatrix of m2'B for 1 particle
   *
   * @param m mass of the particle
   * @param dm2 derivative of the squared particle mass
   * @param type type of particle e.g. fermion/boson
   * @return MatDoub of m2'B
   */
  MatDoub
  calc_m2B(const double &m, const double &dm2, const ParticleType &type);

  /**
   * @brief Calculate the souce term S for 1 particle
   *
   * @param m mass of the particle
   * @param dm2 derivative of the squared particle mass
   * @param dth derivative of theta
   * @param d2th second derivative theta
   * @param type type of particle e.g. fermion/boson
   * @return VecDoub source term
   */
  VecDoub calc_source(const double &m,
                      const double &dm2,
                      const double &dth,
                      const double &d2th,
                      const ParticleType &type);

  /**
   * @brief Check if the boundary have enough decaying modes for the solution to
   * exist.
   *
   * @param MtildeM \f$ M \f$ matrix at the z-negative boundary.
   * @param StildeM \f$ S \f$ vector at the z-negative boundary.
   * @param MtildeP \f$ M \f$ matrix at the z-positive boundary.
   * @param StildeP \f$ S \f$ vector at the z-positive boundary.
   */
  void CheckBoundary(const MatDoub &MtildeM,
                     const VecDoub &StildeM,
                     const MatDoub &MtildeP,
                     const VecDoub &StildeP);

  /**
   * @brief Insert a sub matrix into a block diagonal matrix
   *
   * @param full \f$ reference to the block diagonal matrix \f$
   * @param sub \f$ sub matrix to be inserted \f$
   * @param Stilde \f$ position at which we want to insert the sub matrix \f$
   */
  void InsertBlockDiagonal(MatDoub &full, MatDoub &sub, const size_t position);

  /**
   * @brief Calculate the equations
   *
   * @param z distance to the bubble wall
   * @param Mtilde \f$ A^{-1} \Gamma \f$
   * @param Stilde \f$ A^{-1}\right( \deltaC - m^2' B \left) \f$
   */
  void Equations(const double &z, MatDoub &Mtilde, VecDoub &Stilde);

  /**
   * @brief Solve the transport equation using the shooting method. Runge Kutta
   *
   */
  void ShootingMethod();

  /**
   * @brief Solve the transport equation using the relaxation method.
   * Recommended.
   *
   */
  void RelaxationMethod();

  /**
   * @brief Distributes the points used in the relaxation method
   *
   */
  std::vector<double> MakeDistribution(const double xmax, const size_t npoints);

  /**
   * @brief Solve the transport equations.
   *
   */
  void SolveTransportEquation();

  /**
   * @brief Calculate the baryonic assymetry of the Universe
   *
   */
  void CalculateBAU();

  /**
   * @brief
   *
   * @param size
   * @param Particle
   * @param MuOrU
   */
  void PrintTransportEquation(const int &size,
                              const std::string &Particle,
                              const std::string &MuOrU,
                              const double &multiplier = 1);
};
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT