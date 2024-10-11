#pragma once

#include <BSMPT/utility/NumericalIntegration.h>
#include <cmath>
#include <iostream>
#include <memory>

namespace BSMPT
{

enum K_type
{
  N0,
  Rbar,
  K4,
  D0,
  D1,
  D2,
  Q1,
  Q2,
  Qe1,
  Qe2,
  Q8o1,
  Q8o2,
  Q9o1,
  Q9o2,
};

enum P_type
{
  boson,
  fermion,
  antifermion
};

struct Kinfo
{
  const double Tc, vw;
  double gamw;

  Kinfo(const double T_in, const double vw_in) : Tc(T_in), vw(vw_in)
  {
    gamw = 1 / std::sqrt(1 - vw_in * vw_in);
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

class K4int
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;

public:
  K4int(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  double operator()(const double u);
  ~K4int() {};
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

class Q8o1int1
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;
  double w, pwt;

public:
  double pre;
  Q8o1int1(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  void set_w(const double u_in);
  double operator()(const double y);
  ~Q8o1int1() {};
};

class Q8o1int2
{
private:
  Q8o1int1 integrand;

public:
  Q8o1int2(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : integrand(K_in, s_in, x_in) {};
  double operator()(const double u);
  ~Q8o1int2() {};
};

class Q8o2int1
{
private:
  std::shared_ptr<Kinfo> Ki;
  const double s;
  const double x;
  double w;

public:
  double pre, pwt;
  Q8o2int1(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  void set_w(const double u_in);
  double operator()(const double y);
  ~Q8o2int1() {};
};

class Q8o2int2
{
private:
  Q8o2int1 integrand;

public:
  Q8o2int2(std::shared_ptr<Kinfo> K_in, const double s_in, const double x_in)
      : integrand(K_in, s_in, x_in) {};
  double operator()(const double u);
  ~Q8o2int2() {};
};

class Q9o1int1
{
private:
  std::shared_ptr<Kinfo> Ki;
  const int part;
  const double s;
  const double x;
  double w, pwt;

public:
  double pre2, pre, pre3;
  Q9o1int1(std::shared_ptr<Kinfo> K_in,
           const double s_in,
           const double x_in,
           const int part_in)
      : part(part_in)
      , s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  void set_w(const double u_in);
  double operator()(const double y);
  ~Q9o1int1() {};
};

class Q9o1int2
{
private:
  Q9o1int1 integrand;

public:
  Q9o1int2(std::shared_ptr<Kinfo> K_in,
           const double s_in,
           const double x_in,
           const int part_in)
      : integrand(K_in, s_in, x_in, part_in) {};
  double operator()(const double u);
  ~Q9o1int2() {};
};

class Q9o2int1
{
private:
  std::shared_ptr<Kinfo> Ki;
  const int part;
  const double s;
  const double x;
  double w, pwt;

public:
  double pre;
  Q9o2int1(std::shared_ptr<Kinfo> K_in,
           const double s_in,
           const double x_in,
           const int part_in)
      : part(part_in)
      , s(s_in)
      , x(x_in)
  {
    Ki = K_in;
  }
  void set_w(const double u_in);
  double operator()(const double y);
  ~Q9o2int1() {};
};

class Q9o2int2
{
private:
  Q9o2int1 integrand;

public:
  Q9o2int2(std::shared_ptr<Kinfo> K_in,
           const double s_in,
           const double x_in,
           const int part_in)
      : integrand(K_in, s_in, x_in, part_in) {};
  double operator()(const double u);
  ~Q9o2int2() {};
};
class Kfactor
{
private:
  std::shared_ptr<Kinfo> Ki;
  const bool fast;

public:
  Kfactor(std::shared_ptr<Kinfo> K_in, const bool fast_in) : fast(fast_in)
  {
    Ki = K_in;
  }
  double operator()(const K_type ktype, const P_type ptype, const double m);

  ~Kfactor() {};
};

} // namespace BSMPT
