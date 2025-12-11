// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/baryo_fhck/difeq.h>
#include <BSMPT/utility/utility.h>
#include <fstream>

using BSMPT::Delta;

struct Difeq_DomainWall : Difeq
{
  // Dimension of VEV space
  const size_t &dim;
  // List of z coordinates on the grid
  VecDoub &z;
  // Gradient
  const std::function<std::vector<double>(std::vector<double>)> &dV;
  // Hessian
  const std::function<std::vector<std::vector<double>>(std::vector<double>)>
      &Hessian;

  Difeq_DomainWall(
      // Dimension of VEV space
      const size_t &dim_In,
      // z list
      VecDoub &z_In,
      // Gradient
      const std::function<std::vector<double>(std::vector<double>)> &dV_In,
      // Hessian
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian_In)
      : dim(dim_In)
      , z(z_In)
      , dV(dV_In)
      , Hessian(Hessian_In)
  {
  }

  void smatrix(const int k,
               const int k1,
               const int k2,
               const int jsf,
               VecInt &indexv,
               MatDoub &s,
               MatDoub &y)
  {
    s.zero(); // Set matrix s = 0
    if (k == k1)
    {
      // Boundary conditions dv/dz = 0 on first boundary
      for (size_t field = 0; field < dim; field++)
      {
        // Sn at the first boundary
        s[dim + field][2 * dim + indexv[field]] = 1.0;
        // B0
        s[dim + field][jsf] = y[field][0];
      }
    }
    else if (k > k2 - 1)
    {
      // Boundary conditionsd dv/dz = 0 on second boundary
      for (size_t field = 0; field < dim; field++)
      {
        // Sn at the last boundary
        s[field][2 * dim + indexv[field]] = 1.0;
        // C0
        s[field][jsf] = y[field][k2 - 1];
      }
    }
    else
    {
      const double dz = z[k] - z[k - 1];
      std::vector<double> point; // middle point (phi(k) + phi(k-1)) / 2
      for (size_t j = 0; j < dim; j++)
        point.push_back((y[dim + j][k] + y[dim + j][k - 1]) / 2);
      const std::vector<double> dv                   = dV(point);
      const std::vector<std::vector<double>> hessian = Hessian(point);

      // 0 <= j < dim -> phi'(z)
      for (size_t j = 0; j < dim; j++)
      {
        for (size_t n = 0; n < dim; n++)
        {
          s[j][0 * dim + indexv[n]] = -Delta(j, n);                  // 1
          s[j][1 * dim + indexv[n]] = -1. / 2. * dz * hessian[j][n]; // 3
          s[j][2 * dim + indexv[n]] = Delta(j, n);                   // 5
          s[j][3 * dim + indexv[n]] = -1. / 2. * dz * hessian[j][n]; // 7
        }
        //  Equations for E(k,k-1)
        s[j][jsf] = y[j][k] - y[j][k - 1] - dz * dv[j];
      }

      // dim <= j < 2 * dim -> phi(z)
      for (size_t j = 0; j < dim; j++)
      {
        for (size_t n = 0; n < dim; n++)
        {
          s[dim + j][0 * dim + indexv[n]] = -dz / 2. * Delta(j, n); // 2
          s[dim + j][1 * dim + indexv[n]] = -Delta(j, n);           // 4
          s[dim + j][2 * dim + indexv[n]] = -dz / 2. * Delta(j, n); // 6
          s[dim + j][3 * dim + indexv[n]] = Delta(j, n);            // 8
        }
        //  Equations for E(k,k-1)
        s[dim + j][jsf] = y[dim + j][k] - y[dim + j][k - 1] -
                          dz * (y[j][k] + y[j][k - 1]) / 2;
      }
    }
  }
};