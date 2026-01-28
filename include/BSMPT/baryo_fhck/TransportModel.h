#pragma once

#include <BSMPT/Kfactors/Kernels.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
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

class TransportModel
{
private:
  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief EtaInterface object. To  use BSMPTv2 functions
   *
   */
  std::shared_ptr<CalculateEtaInterface> EtaInterface;

  /**
   * @brief Transition temperature
   *
   */
  double Tstar;

  /**
   * @brief pointer to the vacuum profile solver
   *
   */
  std::shared_ptr<VacuumProfileNS::VacuumProfile> vacuumprofile;

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
   * @brief Coex phase object
   *
   */
  std::shared_ptr<CoexPhases> CoexPhase;

  /**
   * @brief Splines to store the variation of the quark masses
   *
   */
  std::vector<tk::spline> QuarkMassesRe, QuarkMassesIm;

public:
  /**
   * @brief The standard 4 particle that appear in the transport equations
   *
   */
  std::vector<ParticleType> involvedparticles = {ParticleType::LeftFermion,
                                                 ParticleType::RightFermion,
                                                 ParticleType::LeftFermion,
                                                 ParticleType::Boson};
  /**
   * @brief Wall thickness
   *
   */
  double Lw;

  /**
   * @brief Bubble wall velocity
   *
   */
  double vwall;

  /**
   * @brief Mode of the vev profile.
   * - Kink -> Kink solution
   * - FieldEquation -> Calculate the vev profile
   * using the tunneling solution
   */
  VevProfileMode VevProfile = VevProfileMode::Unset;

  /**
   * @brief Construct a new Transport Model object
   *
   * @param pointer_in  Model pointer
   * @param FalseVacuum_In False Vacuum
   * @param TrueVacuum_In True Vacum
   * @param vwall_in Wall velocity
   * @param Tstar_in  Transition temperature
   * @param VevProfile_In Solver mode. Default: field EOM solution
   */
  TransportModel(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const std::vector<double> FalseVacuum_In,
      const std::vector<double> TrueVacuum_In,
      const double &vwall_in,
      const double &Tstar_in,
      const VevProfileMode &VevProfile_In = VevProfileMode::FieldEquation);

  /**
   * @brief Construct a new Transport Model object
   *
   * @param pointer_in Model pointer
   * @param CoexPhase_in Coexphase pointer
   * @param vwall_in Wall velocity
   * @param Tstar_in Transition temperature
   * @param VevProfile_In Solver mode.
   */
  TransportModel(const std::shared_ptr<Class_Potential_Origin> &pointer_in,
                 const std::shared_ptr<CoexPhases> &CoexPhase_in,
                 const double &vwall_in,
                 const double &Tstar_in,
                 const VevProfileMode &VevProfile_In);

  /**
   * @brief Create the VEV vectors
   *
   */
  void Initialize();

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
   * @param zList z distribution where we interpolate the fermion mass
   */
  void GenerateFermionMass(const std::vector<double> &zList);

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
  double GetWMass(const double &z, const double &T);

  /**
   * @brief This function calculates the EW breaking VEV from all contributing
   * field configurations.
   *
   * @param z distance from the bubble wall
   */
  double EWSBVEV(const double &z);

  ~TransportModel() {};
};

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT