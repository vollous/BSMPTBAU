#pragma once

#include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/baryo_fhck/solvde.h>
#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/bounce_solution/bounce_solution.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/utility/utility.h>
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
  TunnelPath
};

class TransportEquations
{
public:
  /**
   * @brief \f$ \eta = \frac{n_B}{n_\gamma}\f$
   *
   */
  std::optional<double> BAUEta;

  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief Transition temperature
   *
   */
  double Tstar;

  /**
   * @brief Bubble wall velocity
   *
   */
  double vwall;

  /**
   * @brief Wall thickness
   *
   */
  std::optional<double> Lw;

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
   * @brief Store the Zs from the solution from the relaxation method
   *
   */
  std::optional<VecDoub> SolutionZ;
  /**
   * @brief Store the solution from the relaxation method
   *
   */
  std::optional<MatDoub> Solution;

  /**
   * @brief Number of Weil fermions
   *
   */
  int nFermions = 3;

  /**
   * @brief Number of bosons
   *
   */
  int nBosons = 1;

  /**
   * @brief Number of equations
   *
   */
  size_t nFB2;

  /**
   * @brief Number of steps in space
   *
   */
  int NumberOfSteps = 101;

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
   * - TunnelPath -> Calculate the vev profile using the tunneling solution
   */
  VevProfileMode VevProfile = VevProfileMode::Unset;

  /**
   * @brief Ki object that has the bubble wall properties
   *
   */
  std::shared_ptr<Kinfo> Ki;

  /**
   * @brief Kfac object to calculate the K-functions
   *
   */
  std::shared_ptr<Kfactor> Kfac;

  /**
   * @brief Coex phase object
   *
   */
  std::shared_ptr<CoexPhases> CoexPhase;

  /**
   * @brief Pointer to Action solution from GW routines
   *
   */
  std::shared_ptr<BounceActionInt> ActionInt;

  /**
   * @brief Pointer to Bounce solution from GW routines
   *
   */
  std::shared_ptr<BounceSolution> Bounce;

  /**
   * @brief Gravitational Wave object
   *
   */
  std::shared_ptr<GravitationalWave> GW;

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param pointer_in Model pointer
   * @param vwall_in Wall velocity
   * @param Tstar_in Transition temperature
   * @param FalseVacuum_in False vacuum
   * @param TrueVacuum_in True vacuum
   */
  TransportEquations(const std::shared_ptr<Class_Potential_Origin> &pointer_in,
                     const double &vwall_in,
                     const double &Tstar_in,
                     const std::vector<double> &FalseVacuum_in,
                     const std::vector<double> &TrueVacuum_in);

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param pointer_in Model pointer
   * @param CoexPhase_in Coexphase pointer
   * @param vwall_in Wall velocity
   * @param Tstar_in Transition temperature
   */
  TransportEquations(const std::shared_ptr<Class_Potential_Origin> &pointer_in,
                     const std::shared_ptr<CoexPhases> &CoexPhase_in,
                     const double &vwall_in,
                     const double &Tstar_in);

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param pointer_in Model pointer
   * @param vwall_in Wall velocity
   * @param Tstar_in Transition temperature
   * @param ActionInt_in ActionInt pointer
   * @param Mode_in Mode. Either Kink of TunnelPath
   */
  TransportEquations(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const double &vwall_in,
      const double &Tstar_in,
      const std::shared_ptr<BounceActionInt> &ActionInt_in,
      const VevProfileMode &Mode_in = VevProfileMode::TunnelPath);

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param pointer_in Model pointer
   * @param Bounce_in BounceSolution pointer
   * @param which_transition_temp Which transition temperature to use ()
   * @param Mode_in Mode. Either Kink of TunnelPath
   */
  TransportEquations(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const std::shared_ptr<BounceSolution> &Bounce_in,
      const int &which_transition_temp,
      const VevProfileMode &Mode_in = VevProfileMode::TunnelPath);

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param pointer_in Model pointer
   * @param GW_in GravitationalWave pointer
   * @param Mode_in Mode. Either Kink of TunnelPath
   */
  TransportEquations(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const std::shared_ptr<GravitationalWave> &GW_in,
      const VevProfileMode &Mode_in = VevProfileMode::TunnelPath);

  /**
   * @brief Create the VEV vectors and initalize **Ki()** and **Kfac()**
   * objects.
   *
   */
  void Initialize();

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
   * @brief Get the Fermion Mass object Calculate the fermion mass and its
   * derivatives
   *
   * @param z distance to the bubble wall
   * @param fermion which fermion, 0 = most massive
   * @param esQuark Quark mass matrix
   * @param m2 \f$ m^2 \f$
   * @param m2prime \f$ m'^2 \f$
   * @param thetaprime \f$ \theta' \f$
   * @param theta2prime \f$ \theta'' \f$
   */
  void
  GetFermionMass(const double &z,
                 const int &fermion,
                 const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> &esQuark,
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
  MatDoub calc_Ainv(const double &m, const P_type &type);

  /**
   * @brief Calculate the 2x2 submatrix of m2'B for 1 particle
   *
   * @param m mass of the particle
   * @param dm2 derivative of the squared particle mass
   * @param type type of particle e.g. fermion/boson
   * @return MatDoub of m2'B
   */
  MatDoub calc_m2B(const double &m, const double &dm2, const P_type &type);

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
                      const P_type &type);

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
  VecDoub MakeDistribution(const double xmax, const size_t npoints);

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
                              const double &multiplier = 3);
};
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT