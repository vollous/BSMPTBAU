#pragma once

#include <BSMPT/utility/NumericalIntegration.h>
#include <cmath>
#include <iostream>
#include <memory>

namespace BSMPT
{

enum K_type
{
  N0bos,
  N0fer,
  Rbarbos,
  Rbarfer,
  D0bos,
  D0fer,
  D1bos,
  D1fer,
  D2bos,
  D2fer,
  Q1bos,
  Q1fer,
  Q2bos,
  Q2fer,
  Qe1bos,
  Qe1fer,
  Qe2bos,
  Qe2fer,
  Q8o1bos,
  Q8o1fer,
  Q8o2bos,
  Q8o2fer,
  Q9o1bos,
  Q9o1fer,
  Q9o2bos,
  Q9o2fer
};

struct Kinfo
{
  const double Tc, vw;
  double gamw;
  int n, m, k;
  const bool fast;

  Kinfo(const double T_in, const double vw_in, const bool fast_in)
      : Tc(T_in)
      , vw(vw_in)
      , fast(fast_in)
  {
    gamw = 1 / std::sqrt(1 - vw_in * vw_in);
  }

  void set_nmk(const double n_in, const double m_in, const double k_in)
  {
    n = n_in;
    m = m_in;
    k = k_in;
  }

  double Kintegrand2D(const double w, const double y, const double x)
  {
    double pwt = sqrt(w * w - x * x);
    double pzt = gamw * (y * pwt - w * vw);
    double Et  = gamw * (w - vw * y * pwt);
    double Vx  = pow((pzt / sqrt(pzt * pzt + x * x)), 2) * 1 /
                sqrt(1 - pzt * pzt / (Et * Et));
    return -3 / (M_PI * M_PI * gamw) * pow(Tc, n - m - k + 1) * Vx * pwt *
           pow(pzt, n) / pow(Et, m - 1);
  }
  ~Kinfo() {};
};

double f0w(const double w, const double s, const int diff);

class N0int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  N0int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~N0int() {};
};

class Rbarint
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  Rbarint(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~Rbarint() {};
};

class D0int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  D0int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~D0int() {};
};

class D2int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  D2int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~D2int() {};
};

class Q1int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  Q1int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~Q1int() {};
};

class Q2int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  Q2int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~Q2int() {};
};

class Qe1int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  Qe1int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~Qe1int() {};
};

class Qe2int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  Qe2int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  };
  double operator()(const double u);
  ~Qe2int() {};
};

class Q8oint1
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;
  double u;

public:
  Q8oint1(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  void set_u(const double u_in);
  double operator()(const double y);
  ~Q8oint1() {};
};

class Q8oint2
{
private:
  Q8oint1 integrand;

public:
  Q8oint2(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : integrand(K_in, s_in, x_in) {};
  double operator()(const double u);
  ~Q8oint2() {};
};

class Q9oint1
{
private:
  std::shared_ptr<Kinfo> Ki;
  const int part;
  const double s;
  const double x;
  double u;

public:
  Q9oint1(std::shared_ptr<Kinfo> K_in,
          const double s_in,
          const double x_in,
          const int part_in)
      : part(part_in)
      , s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  void set_u(const double u_in);
  double operator()(const double y);
  ~Q9oint1() {};
};

class Q9oint2
{
private:
  Q9oint1 integrand;

public:
  Q9oint2(std::shared_ptr<Kinfo> K_in,
          const double s_in,
          const double x_in,
          const int part_in)
      : integrand(K_in, s_in, x_in, part_in) {};
  double operator()(const double u);
  ~Q9oint2() {};
};

class Kfactor
{
private:
  std::shared_ptr<Kinfo> Ki;

public:
  Kfactor(std::shared_ptr<Kinfo> K_in) { Ki = K_in; }

  double operator()(const K_type type, const double m);

  ~Kfactor() {};
};

} // namespace BSMPT
