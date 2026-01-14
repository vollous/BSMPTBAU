#pragma once

#include <BSMPT/utility/NumericalIntegration.h>
#include <cmath>
#include <iostream>
#include <memory>

namespace BSMPT
{

enum KernelType
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

enum ParticleType
{
  Boson,
  Fermion,
};

double f0(const double x, const double s, const int diff);

class KernelInty
{
private:
  const int k, structure, n, m;
  const double vw, gamw, s, x;
  double w, pwt;

public:
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
  void set_all(const double u);
  double operator()(const double y);
  ~KernelInty() {};
};

class KernelIntw
{
private:
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
  double operator()(const double u);
  ~KernelIntw() {};
};

class Q9KernelInty
{
private:
  const int structure, l;
  const double vw, gamw, s, x;
  double w, pwt;

public:
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
  void set_all(const double u);
  double operator()(const double y);
  ~Q9KernelInty() {};
};

class Q9KernelIntw
{
private:
  Q9KernelInty Integrand;

public:
  Q9KernelIntw(const int structure_in,
               const int l_in,
               const double vw_in,
               const double gamw_in,
               const double s_in,
               const double x_in)
      : Integrand(structure_in, l_in, vw_in, gamw_in, s_in, x_in) {};
  double operator()(const double u);
  ~Q9KernelIntw() {};
};

class N0Int
{
private:
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
  double operator()(const double u);
  ~N0Int() {};
};

class RbarInt
{
private:
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
  double operator()(const double u);
  ~RbarInt() {};
};

class K4Int
{
private:
  const double s, x;

public:
  K4Int(const double s_in, const double x_in) : s(s_in), x(x_in) {}
  double operator()(const double u);
  ~K4Int() {};
};

class Kernel
{
private:
  const int l;         // lth - moment
  const int structure; // 1 | V = Vs, 2 | V = Vh

public:
  Kernel(const int l_in, const int structure_in)
      : l(l_in)
      , structure(structure_in) {};
  double operator()(const KernelType Kern,
                    const ParticleType P,
                    const double x,
                    const double vw);
  ~Kernel() {};
};

} // namespace BSMPT
