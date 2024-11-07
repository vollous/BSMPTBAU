#pragma once

#include <BSMPT/utility/matrix_operations.h>
#include <cassert>
#include <cmath>
#include <iostream>

class LUdcmp
{
private:
  int n;
  MatDoub lu;
  VecInt indx;

public:
  LUdcmp(MatDoub &a);
  void solve(VecDoub &b, VecDoub &x);
  ~LUdcmp() {};
};

template <class F>
void lnsrch(VecDoub &xold,
            const double fold,
            VecDoub &g,
            VecDoub &p,
            VecDoub &x,
            double &f,
            const double stpmax,
            bool &check,
            F &func)
{
  static const double ALF = 1e-4, TOLX = 1e-10;
  double a, alam, alam2 = 0., alamin, b, disc, f2 = 0.;
  double rhs1, rhs2, slope = 0., sum = 0., temp, test, tmplam;
  int i, n = xold.size();
  check = false;
  for (i = 0; i < n; i++)
    sum += p[i] * p[i];
  sum = std::sqrt(sum);
  if (sum > stpmax)
    for (i = 0; i < n; i++)
      p[i] *= stpmax / sum;
  for (i = 0; i < n; i++)
    slope += g[i] * p[i];
  assert(slope < 0.);
  test = 0.; // this definition can be pulled up
  for (i = 0; i < n; i++)
  {
    temp = std::abs(p[i]) / std::max(std::abs(xold[i]), 1.);
    if (temp > test) test = temp; // not sure why we test this. Can remove?
  }
  alamin = TOLX / test;
  alam   = 1.;
  while (true)
  {
    for (i = 0; i < n; i++)
      x[i] = xold[i] + alam * p[i];
    f = func(x);
    if (alam < alamin)
    {
      for (i = 0; i < n; i++)
        x[i] = xold[i];
      check = true;
      return;
    }
    else if (f <= fold + ALF * alam * slope)
      return;
    else
    {
      if (alam == 1.)
        tmplam = -slope / (2. * (f - fold - slope));
      else
      {
        rhs1 = f - fold - alam * slope;
        rhs2 = f2 - fold - alam2 * slope;
        a    = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
        b    = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) /
            (alam - alam2);
        if (a == 0.)
          tmplam = -slope / (2. * b);
        else
        {
          disc = b * b - 3. * a * slope;
          if (disc < 0.)
            tmplam = 0.5 * alam;
          else if (b <= 0.)
            tmplam = (-b + std::sqrt(disc)) / (3. * a);
          else
            tmplam = -slope / (b + std::sqrt(disc));
        }
        if (tmplam > 0.5 * alam) tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2    = f;
    alam  = std::max(tmplam, 0.1 * alam);
  }
}

template <class F> class NRfdjac
{
private:
  F &func;
  const double EPS = 1e-5;

public:
  NRfdjac(F &funcc) : func(funcc) {};
  MatDoub operator()(VecDoub &x, VecDoub &fvec)
  {
    int n = x.size();
    MatDoub df(n, VecDoub(n));
    VecDoub xh = x;
    for (int j = 0; j < n; j++)
    {
      double temp = xh[j];
      double h    = EPS * std::abs(temp);
      if (h == 0.) h = EPS;
      xh[j]     = temp + h;
      h         = xh[j] - temp;
      VecDoub f = func(xh);
      xh[j]     = temp;
      for (int i = 0; i < n; i++)
        df[i][j] = (f[i] - fvec[i]) / h;
    }
    return df;
  }
  ~NRfdjac() {};
};

template <class F> class NRfmin
{
private:
  F &func;

public:
  VecDoub fvec;
  NRfmin(F &funcc) : func(funcc) {};
  double operator()(VecDoub &x)
  {
    int n      = x.size();
    double sum = 0.;
    fvec       = func(x);
    for (int i = 0; i < n; i++)
      sum += (fvec[i] * fvec[i]);
    return 0.5 * sum;
  }
  ~NRfmin() {};
};

template <class F> void newt(VecDoub &x, bool &check, F &vecfunc)
{
  const int MAXITS  = 200;
  const double TOLF = 1e-5, TOLMIN = 1e-8, STPMAX = 100.;
  const double TOLX = 1e-10;
  int i, j, its, n = x.size();
  double den, f, fold, stpmax, sum, temp, test;
  VecDoub g(n), p(n), xold(n);
  MatDoub fjac(n, VecDoub(n));
  NRfmin<F> fmin(vecfunc);
  NRfdjac<F> fdjac(vecfunc);
  VecDoub &fvec = fmin.fvec;
  f             = fmin(x);
  test          = 0.;
  for (i = 0; i < n; i++)
    if (std::abs(fvec[i]) > test) test = std::abs(fvec[i]);
  if (test < 0.01 * TOLF)
  {
    check = false;
    return;
  }
  sum = 0.;
  for (i = 0; i < n; i++)
    sum += x[i] * x[i];
  stpmax = STPMAX * std::max(std::sqrt(sum), (double)n);
  for (its = 0; its < MAXITS; its++)
  {
    fjac = fdjac(x, fvec);
    for (i = 0; i < n; i++)
    {
      sum = 0.0;
      for (j = 0; j < n; j++)
        sum += fjac[j][i] * fvec[j];
      g[i] = sum;
    }
    for (i = 0; i < n; i++)
      xold[i] = x[i];
    fold = f;
    for (i = 0; i < n; i++)
      p[i] = -fvec[i];
    LUdcmp alu(fjac);
    alu.solve(p, p);
    lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
    test = 0.0;
    for (i = 0; i < n; i++)
      if (std::abs(fvec[i]) > test) test = std::abs(fvec[i]);
    if (test < TOLF)
    {
      check = false;
      return;
    }
    if (check)
    {
      test = 0.0;
      den  = std::max(f, 0.5 * n);
      for (i = 0; i < n; i++)
      {
        temp = std::abs(g[i]) * std::max(std::abs(x[i]), 1.0) / den;
        if (temp > test) test = temp;
      }
      check = (test < TOLMIN);
      return;
    }
    test = 0.0;
    for (i = 0; i < n; i++)
    {
      temp = (std::abs(x[i] - xold[i])) / std::max(std::abs(x[i]), 1.0);
      if (temp > test) test = temp;
    }
    if (test < TOLX) return;
  }
}
