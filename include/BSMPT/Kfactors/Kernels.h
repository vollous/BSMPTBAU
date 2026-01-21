#pragma once

#include <BSMPT/utility/NumericalIntegration.h>
#include <cmath>
#include <iostream>
#include <memory>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

/**
 * @brief Types of Kernel functions
 *
 */
enum class KernelType
{
  K,
  D,
  Q,
  Qe,
  Q8o,
  Q9o,
  Rb,
  K4FH
};

/**
 * @brief Types of particles
 *
 */
enum class ParticleType
{
  Boson,
  Fermion,
};

/**
 * @brief Function for integer powers. Better suitable than the std one
 *
 * @param x Value of the number
 * @param i Power of the number
 */
double ipow(const double x, const int n);

/**
 * @brief Fermi/Dirac and Bose/Einstein distribution with their derivatives
 *
 * @param x Value at which we want to evaluate
 * @param s Spin statistic. 1 for fermions, -1 for bosons
 * @param diff Order of derivative {0,1,2}
 */
double f0(const double x, const double s, const int diff);

class KernelInty
{
private:
  /**
   * @brief Settings which adjust the kernel function F(f^(k), V, n, m)
   * see [2510.21915]. Structure = {0, 1, 2} -> V={1, V_s, V_h}
   *
   */
  const int k, structure, n, m;

  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

  /**
   * @brief Variables which occur during the integration over y
   *
   */
  double w, pwt;

public:
  /**
   * @brief Useful prefactor to save computation time
   *
   */
  double pre;

  KernelInty(const int k_in,
             const int structure_in,
             const int n_in,
             const int m_in,
             const double vw_in,
             const double gamw_in,
             const double s_in,
             const double x_in)
      : k(k_in)
      , structure(structure_in)
      , n(n_in)
      , m(m_in)
      , vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in) {};

  /**
   * @brief Sets all relevant variables of the class
   * @param  u Integration variable. Defined through w = x + (1 - u) / u. Takes
   * value between 0 and 1
   */
  void set_all(const double u);

  /**
   * @brief operator which evaluates the integrand of the y-dependent part of
   * the kernel integrand
   * @param y Integration variable between [-1, 1]
   */
  double operator()(const double y);
  ~KernelInty() {};
};

class KernelIntw
{
private:
  /**
   * @brief Integrand depending on y.
   */
  KernelInty Integrand;

public:
  KernelIntw(const int k_in,
             const int structure_in,
             const int n_in,
             const int m_in,
             const double vw_in,
             const double gamw_in,
             const double s_in,
             const double x_in)
      : Integrand(k_in, structure_in, n_in, m_in, vw_in, gamw_in, s_in, x_in) {
      };

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);

  ~KernelIntw() {};
};

/**
 * @brief Separate class for Q9 kernels as they involve two normal kernels
 */
class Q9KernelInty
{
private:
  /**
   * @brief Structure = {0, 1, 2} -> V={1, V_s, V_h}. l is the moment of the
   * kernel
   */
  const int structure, l;

  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

  /**
   * @brief Variables which occur during the integration over y
   *
   */
  double w, pwt;

public:
  /**
   * @brief Useful prefactors to save computation time
   *
   */
  double pre1, pre2;

  Q9KernelInty(const int structure_in,
               const int l_in,
               const double vw_in,
               const double gamw_in,
               const double s_in,
               const double x_in)
      : structure(structure_in)
      , l(l_in)
      , vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in) {};

  /**
   * @brief Sets all relevant variables of the class
   * @param  u Integration variable. Defined through w = x + (1 - u) / u. Takes
   * value between 0 and 1
   */
  void set_all(const double u);

  /**
   * @brief operator which evaluates the integrand of the y-dependent part of
   * the kernel integrand
   * @param y Integration variable between [-1, 1]
   */
  double operator()(const double y);

  ~Q9KernelInty() {};
};

class Q9KernelIntw
{
private:
  /**
   * @brief Integrand depending on y.
   */
  Q9KernelInty Integrand;

public:
  Q9KernelIntw(const int structure_in,
               const int l_in,
               const double vw_in,
               const double gamw_in,
               const double s_in,
               const double x_in)
      : Integrand(structure_in, l_in, vw_in, gamw_in, s_in, x_in) {};

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);
  ~Q9KernelIntw() {};
};

class N0Int
{
private:
  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

public:
  N0Int(const double vw_in,
        const double gamw_in,
        const double s_in,
        const double x_in)
      : vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in)
  {
  }

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);
  ~N0Int() {};
};

class RbarInt
{
private:
  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

public:
  RbarInt(const double vw_in,
          const double gamw_in,
          const double s_in,
          const double x_in)
      : vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in) {};

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);
  ~RbarInt() {};
};

class K4Int
{
private:
  /**
   * @brief spin-statistic, lower boundary of the integral
   *
   */
  const double s, x;

public:
  K4Int(const double s_in, const double x_in) : s(s_in), x(x_in) {}

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);

  ~K4Int() {};
};

class Kernel
{
private:
  /**
   * @brief Moment of the Kernel
   */
  const int l;

  /**
   * @brief Structure = {0, 1, 2} -> V={1, V_s, V_h}
   */
  const int structure;

public:
  Kernel(const int l_in, const int structure_in)
      : l(l_in)
      , structure(structure_in) {};
  /**
   * @brief Evaluates the needed kernel at a given x and wall velocity
   * @param KernelType Type of Kernel
   * @param ParticleType Type of Particle
   * @param x Value at which we want to evaluate the kernel x=m/T
   * @param vw Wall velocity at which we want to evaluate the kernel
   */
  double operator()(const KernelType Kern,
                    const ParticleType P,
                    const double x,
                    const double vw);
  ~Kernel() {};
};

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT