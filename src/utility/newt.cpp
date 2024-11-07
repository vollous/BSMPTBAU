#include <BSMPT/utility/newt.h>

LUdcmp::LUdcmp(MatDoub &a) : n(a.size()), lu(a), indx(n)
{
  static const double TINY = 1.0e-40;
  int i, imax, j, k;
  double big, temp;
  VecDoub vv(n);
  for (i = 0; i < n; i++)
  {
    big = 0.;
    for (j = 0; j < n; j++)
    {
      temp = std::abs(lu[i][j]);
      if (temp > big) big = temp;
    }
    assert(big != 0.);
    vv[i] = 1. / big;
  }
  for (k = 0; k < n; k++)
  {
    big = 0.;
    for (i = k; i < n; i++)
    {
      temp = vv[i] * std::abs(lu[i][k]);
      if (temp > big)
      {
        big  = temp;
        imax = i;
      }
    }
    if (k != imax)
    {
      for (j = 0; j < n; j++)
      {
        temp        = lu[imax][j];
        lu[imax][j] = lu[k][j];
        lu[k][j]    = temp;
      }
      vv[imax] = vv[k];
    }
    indx[k] = imax;
    if (lu[k][k] == 0.) lu[k][k] = TINY;
    for (i = k + 1; i < n; i++)
    {
      temp = lu[i][k] /= lu[k][k];
      for (j = k + 1; j < n; j++)
        lu[i][j] -= temp * lu[k][j];
    }
  }
}

void LUdcmp::solve(VecDoub &b, VecDoub &x)
{
  int i, ii = 0, ip, j;
  double sum;
  assert(b.size() == n && x.size() == n);
  for (i = 0; i < n; i++)
    x[i] = b[i];
  for (i = 0; i < n; i++)
  {
    ip    = indx[i];
    sum   = x[ip];
    x[ip] = x[i];
    if (ii != 0)
      for (j = ii - 1; j < i; j++)
        sum -= lu[i][j] * x[j];
    else if (sum != 0.)
      ii = i + 1;
    x[i] = sum;
  }
  for (i = n - 1; i >= 0; i--)
  {
    sum = x[i];
    for (j = i + 1; j < n; j++)
      sum -= lu[i][j] * x[j];
    x[i] = sum / lu[i][i];
  }
}