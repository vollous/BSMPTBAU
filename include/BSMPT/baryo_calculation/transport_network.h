#pragma once

#include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/NumericalIntegration.h>
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
  std::shared_ptr<Kinfo> Ki;
  Kfactor K;
  double LW, vh, vs, LAM = 1000.;
  double D0b, D1b, D2b, D1h, D2h;
  double detAb, detAh;
  double Gtotb, Gtoth;

public:
  TransportNetwork(std::shared_ptr<Kinfo> K_in) : K(K_in)
  {
    Ki    = K_in;
    LW    = 5. / Ki->Tc;
    vh    = Ki->Tc;
    vs    = 2 * Ki->Tc;
    D0b   = K(D0, fermion, 0.);
    D1b   = K(D1, fermion, 0.);
    D2b   = K(D2, fermion, 0.);
    D1h   = K(D1, boson, 0.);
    D2h   = K(D2, boson, 0.);
    detAb = Ki->vw * D1b + D2b;
    detAh = Ki->vw * D1h + D2h;
    Gtotb = K(K4, fermion, 0.) * Ki->Tc / (6. * D0b);
    Gtoth = K(K4, boson, 0.) * Ki->Tc / (20. * K(D0, boson, 0.));
  }

  P_type get_particle_type(const Particles prtcl);

  double vev1_profile(const double &z, const size_t deriv);

  double vev2_profile(const double &z, const size_t deriv);

  double top_mass(const double z, const size_t deriv);

  double W_mass(const double z);

  double theta(const double z, const size_t deriv);

  MatDoub Ainv(const double z);

  MatDoub gamma(const double z);

  void
  spline_Kfactors(const double zmin, const double zmax, const size_t N_points);

  VecDoub calc_Source(const double z);

  void operator()(const double z, VecDoub &u, VecDoub &du);

  ~TransportNetwork() {};
};

class dgl
{
private:
  /* data */
public:
  dgl(/* args */) {};
  void operator()(const double t, const VecDoub &y, VecDoub &dy);
  ~dgl() {};
};

class shootf
{
private:
  // const double zl = -0.7, zr = 0.9;
  TransportNetwork tr;
  dgl eq;
  const double xi = 0., xf = 10.;
  const double Nb  = 1.;
  const double Neq = 2.;
  const double Np  = 20.;
  double dx;

public:
  bool save   = false;
  double zmid = 0.;
  shootf(std::shared_ptr<Kinfo> K_in) : tr(K_in)
  {
    dx = (xf - xi) / Np;
    // tr.spline_Kfactors(zl, zr, 600);
    // MatDoub a = tr.Ainv(zr);
    // MatDoub b = tr.gamma(zr);
    // printmat(a * b);
    // exit(1);
  };
  VecDoub operator()(VecDoub &v);
  ~shootf() {};
};

class shootf2
{
private:
  const double zl = -0.7, zr = 0.8;
  TransportNetwork tr;

public:
  bool save   = false;
  double zmid = 0.;
  shootf2(std::shared_ptr<Kinfo> K_in) : tr(K_in)
  {
    tr.spline_Kfactors(zl, zr, 600);
    // MatDoub a = tr.Ainv(zr);
    // MatDoub b = tr.gamma(zr);
    // printmat(a * b);
    // exit(1);
  };
  VecDoub operator()(VecDoub &v);
  ~shootf2() {};
};

} // namespace BSMPT
