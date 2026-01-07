#pragma once

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/SMparam.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/difeq_vacuum_profile.h>
#include <iostream>
#include <string>
#include <vector>

namespace BSMPT
{
namespace VacuumProfileNS
{

struct VacuumProfile
{
  /**
   * @brief Method used to solve the field equation
   * - Field - Fixed \f$ \vec\phi(z \to \infty) = \vec\phi_f \f$ and \f$
   * \vec\phi(z \to -\infty) = \vec\phi_t \f$
   * - Deriv - Fixes \frac{d\vec\phi}{dz}(|z| \to \infty) = 0
   *
   */
  ProfileSolverMode mode = ProfileSolverMode::Deriv;
  /**
   * @brief Dimension of VEV space
   *
   */
  const size_t &dim;
  /**
   * @brief True and False Vacuum
   *
   */
  const std::vector<double> &TrueVacuum, &FalseVacuum;
  /**
   * @brief Potential
   *
   */
  const std::function<double(std::vector<double>)> &V;
  /**
   * @brief Gradient
   *
   */
  const std::function<std::vector<double>(std::vector<double>)> &dV;
  /**
   * @brief Hessian
   *
   */
  const std::function<std::vector<std::vector<double>>(std::vector<double>)>
      &Hessian;
  /**
   * @brief z positions corresponding to the path knots
   *
   */
  std::vector<double> z;
  /**
   * @brief Tunneling
   *
   */
  std::vector<std::vector<double>> &path;

  /**
   * @brief Relocate the vacuum profile "center" to \f$ z = 0 \f$ but putting
   * \f$ \max \left\{\frac{d\vec\phi}{dz}\right\}^2 \f$ at \f$ z = 0 \f$
   *
   */
  void CenterPath();

  /**
   * @brief Calculate the vacuum profile
   *
   */
  void CalculateProfile();

  /**
   * @brief Construct a new Vacuum Profile object. Attemps to calculate \f$ L_w
   * \f$ for the first guess using the kink solution
   *
   * @param dim_In VEV dimension
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   * @param V_In Potential
   * @param dV_In Potential Gradient
   * @param Hessian_In Potential Hessian
   */
  VacuumProfile(
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
          &Hessian_In);

  /**
   * @brief Construct a new Vacuum Profile object. Uses \f$ L_w \f$ for the
   * first path guess using the kink solution
   *
   * @param dim_In VEV dimension
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   * @param V_In Potential
   * @param dV_In Potential Gradient
   * @param Hessian_In Potential Hessian
   */
  VacuumProfile(
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
      const double &Lw_In);

  /**
   * @brief Construct a new Vacuum Profile object
   *
   * @param dim_In VEV dimension
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   * @param V_In Potential
   * @param dV_In Potential Gradient
   * @param Hessian_In Potential Hessian
   * @param z_In List of the \f$ z \f$ position
   * @param path_In List of the field knots
   */
  VacuumProfile(
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
      std::vector<std::vector<double>> &path_In);
};

} // namespace VacuumProfileNS
} // namespace BSMPT
