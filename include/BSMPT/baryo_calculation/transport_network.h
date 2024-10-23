#pragma once

#include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/matrix_operations.h>
#include <memory>
#include <vector>

enum Particles
{
  tL,
  tR,
  bL,
  bR,
  h,
  W
};

namespace BSMPT
{

class TransportNetwork
{
private:
  std::shared_ptr<Class_Potential_Origin> modelPointer;
  std::shared_ptr<Kinfo> Ki;
  const std::vector<Particles> prtcl_list;
  VecDoub vev_critical;
  Kfactor K;
  const double thsym = 0.0603557;
  const double thbrk = 0.00688404;
  const double LW    = 0.0341292;

public:
  TransportNetwork(std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
                   std::shared_ptr<Kinfo> K_in,
                   std::vector<Particles> prtcl_list_in,
                   const VecDoub vev_critical_in)
      : prtcl_list(prtcl_list_in)
      , K(K_in)
  {
    Ki           = K_in;
    modelPointer = modelPointer_input;
    vev_critical = modelPointer->MinimizeOrderVEV(vev_critical_in);
  }

  size_t get_N_particles();

  P_type get_particle_type(const Particles prtcl);

  double vevProfileKink(const double &z, size_t deriv);

  VecDoub calc_vev(const double &z, size_t deriv);

  double calculate_theta(const double &z, size_t diff);

  VecDoub get_top_mass_and_derivative(const VecDoub &vev) const;

  VecDoub get_h_mass_and_derivative(const std::vector<double> &vev);

  double get_W_mass(const VecDoub &vev) const;

  VecDoub get_squared_mass_and_deriv(const double z, const Particles prtcl);

  void spline_Kfactors(const double zmin, const double zmax, const size_t N_points);

  MatDoub calc_A_inv(const double z);

  MatDoub calc_B(const double z);

  MatDoub calc_Collision(const double z);

  VecDoub calc_Source(const double z);

  void operator()(const double z, VecDoub &u, VecDoub &du);

  ~TransportNetwork() {};
};
} // namespace BSMPT
