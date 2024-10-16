#include <BSMPT/utility/matrix_operations.h>

VecDoub operator*(const VecDoub &a, const double b)
{
  VecDoub res(a.size());
  for (size_t i = 0; i < a.size(); i++)
    res[i] = a[i] * b;
  return res;
}

VecDoub operator*(const double a, const VecDoub &b)
{
  VecDoub res(b.size());
  for (size_t i = 0; i < b.size(); i++)
    res[i] = a * b[i];
  return res;
}

VecDoub operator+(const VecDoub &a, const VecDoub &b)
{
  VecDoub res(a.size());
  for (size_t i = 0; i < a.size(); i++)
    res[i] = a[i] + b[i];
  return res;
}

VecDoub operator*(const MatDoub &a, const VecDoub &b)
{
  assert(a[0].size() == b.size());
  VecDoub res(a[0].size(), 0);
  for (size_t i = 0; i < a.size(); i++)
    for (size_t j = 0; j < a.size(); j++)
      res[i] += a[i][j] * b[j];
  return res;
}

MatDoub operator+(const MatDoub &a, const MatDoub &b)
{
  assert(a.size() == b.size());
  assert(a[0].size() == b[0].size());
  MatDoub res(a.size(), VecDoub(a[0].size(), 0.));
  for (size_t i = 0; i < a.size(); i++)
    for (size_t j = 0; j < a.size(); j++)
      res[i][j] = a[i][j] + b[i][j];
  return res;
}

MatDoub operator-(const MatDoub &a, const MatDoub &b)
{
  assert(a.size() == b.size());
  assert(a[0].size() == b[0].size());
  MatDoub res(a.size(), VecDoub(a[0].size(), 0.));
  for (size_t i = 0; i < a.size(); i++)
    for (size_t j = 0; j < a.size(); j++)
      res[i][j] = a[i][j] - b[i][j];
  return res;
}

MatDoub operator*(const MatDoub &a, const MatDoub &b)
{
  assert(a[0].size() == b.size());
  MatDoub res(a.size(), VecDoub(b[0].size(), 0.));
  for (size_t i = 0; i < a.size(); i++)
    for (size_t j = 0; j < b[0].size(); j++)
      for (size_t k = 0; k < b.size(); k++)
        res[i][j] += a[i][k] * b[k][j];
  return res;
}

void set_zero(MatDoub &a)
{
  for (auto &it : a)
    for (auto &jt : it)
      jt = 0;
}

void printvec(const VecDoub &a)
{
  std::cout << "----------------------------------\n";
  for (auto it : a)
    std::cout << it << "\t";
  std::cout << "\n";
}

void printmat(const MatDoub &a)
{
  std::cout << "----------------------------------\n";
  for (auto it : a)
  {
    for (auto jt : it)
    {
      std::cout << jt << " ";
    }
    std::cout << "\n";
  }
}
