// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/baryo_fhck/difeq.h>
#include <BSMPT/utility/utility.h>
#include <fstream>

using BSMPT::Delta;

enum class ProfileSolverMode
{
  Field,
  Deriv
};

struct Difeq_BubbleProfile : Difeq
{
  // Dimension of VEV space
  const size_t &dim;
  // List of z coordinates on the grid
  VecDoub &z;
  // False and true vacuum
  const std::vector<double> &TrueVacuum, &FalseVacuum;
  // Relxation method "method"
  const ProfileSolverMode &mode;
  // Gradient
  const std::function<std::vector<double>(std::vector<double>)> &dV;
  // Hessian
  const std::function<std::vector<std::vector<double>>(std::vector<double>)>
      &Hessian;

  Difeq_BubbleProfile(
      // Dimension of VEV space
      const size_t &dim_In,
      // z list
      VecDoub &z_In,
      // True
      const std::vector<double> &TrueVacuum_In,
      // False
      const std::vector<double> &FalseVacuum_In,
      // Mode
      const ProfileSolverMode &mode_In,
      // Gradient
      const std::function<std::vector<double>(std::vector<double>)> &dV_In,
      // Hessian
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian_In)
      : dim(dim_In)
      , z(z_In)
      , TrueVacuum(TrueVacuum_In)
      , FalseVacuum(FalseVacuum_In)
      , mode(mode_In)
      , dV(dV_In)
      , Hessian(Hessian_In)
  {
    std::cout << "We have \t" << dim_In << " VEVs\n";
    std::ofstream myfile;
    myfile.open("z_list_profile.tsv");
    for (size_t m = 0; m < z_In.size(); m++)
      myfile << z_In[m] << "\t";
    myfile.close();
  }

  /*void printmatrix(const MatDoub &ss)
  {
    size_t ne = 2 * dim;
    for (size_t i = 0; i < ne; i++)
    {
      for (size_t j = 0; j < 2 * ne + 1; j++)
        std::cout << ss[i][j] << "\t";
      std::cout << "\n";
    }
  }*/

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
      if (mode == ProfileSolverMode::Deriv)
      {
        // Boundary conditions dv/dz = 0 on first boundary
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the first boundary
          s[dim + field][2 * dim + 1 + indexv[field]] = 1.0;
          // B0
          s[dim + field][jsf] = y[field][0];
        }
        // Sn at the first boundary
        s[dim + dim][2 * dim + 1 + indexv[dim]] = 1.0;
        // B0
        s[dim + dim][jsf] = y[dim][0] + 1.0;
      }
      else if (mode == ProfileSolverMode::Field)
      {
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the first boundary
          s[dim + field][2 * dim + indexv[field]] = 1.0;
          // B0
          s[dim + field][jsf] = y[dim + field][0] - TrueVacuum[field];
        }
        // Sn at the first boundary
        s[dim + dim][2 * dim + 1 + indexv[0]] = 1.0;
        // B0
        s[dim + dim][jsf] = y[0][0];
        std::cout << "\nPrint matrix\n";
        printmatrix(s);
        std::cout << "\nPrint matrix\n";
      }
    }
    else if (k > k2 - 1)
    {
      if (mode == ProfileSolverMode::Deriv)
      {
        // Boundary conditionsd dv/dz = 0 on second boundary
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the last boundary
          s[field][2 * dim + 1 + indexv[field]] = 1.0;
          // C0
          s[field][jsf] = y[field][k2 - 1];
        }
      }
      else if (mode == ProfileSolverMode::Field)
      {
        // Fix the fields at |z| -> Infinity
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the last boundary
          s[field][2 * dim + 1 + dim + indexv[field]] = 1.0;
          // C0
          s[field][jsf] = y[dim + field][k2 - 1] - FalseVacuum[field];
        }
      }
    }
    else
    {
      const double dz = z[k] - z[k - 1];
      std::vector<double> point; // middle point (phi(k) + phi(k-1)) / 2
      for (size_t j = 0; j < dim; j++)
        point.push_back((y[dim + j][k] + y[dim + j][k - 1]) / 2);
      std::vector<double> dv                   = dV(point);
      std::vector<std::vector<double>> hessian = Hessian(point);

      // 0 <= j < dim -> phi'(z)
      for (size_t j = 0; j < dim; j++)
      {
        for (size_t n = 0; n < dim; n++)
        {
          s[j][0 * dim + indexv[n]]     = -Delta(j, n);                  // 1
          s[j][1 * dim + indexv[n]]     = -1. / 2. * dz * hessian[j][n]; // 3
          s[j][2 * dim + indexv[n] + 1] = Delta(j, n);                   // 5
          s[j][3 * dim + indexv[n] + 1] = -1. / 2. * dz * hessian[j][n]; // 7
        }
        s[j][1 * dim + indexv[dim]]     = dz / 2. * y[j][k - 1]; // 11
        s[j][3 * dim + 1 + indexv[dim]] = dz / 2. * y[j][k];     // 16
        //  Equations for E(k,k-1)
        s[j][jsf] = y[j][k] - y[j][k - 1] - dz * dv[j];
      }

      // dim <= j < 2 * dim -> phi(z)
      for (size_t j = 0; j < dim; j++)
      {
        for (size_t n = 0; n < dim; n++)
        {
          s[dim + j][0 * dim + indexv[n]]     = -dz / 2. * Delta(j, n); // 2
          s[dim + j][1 * dim + indexv[n]]     = -Delta(j, n);           // 4
          s[dim + j][2 * dim + indexv[n] + 1] = -dz / 2. * Delta(j, n); // 6
          s[dim + j][3 * dim + indexv[n] + 1] = Delta(j, n);            // 8
        }
        s[dim + j][2 * dim]     = 0; // 12
        s[dim + j][4 * dim + 1] = 0; // 17
        //  Equations for E(k,k-1)
        s[dim + j][jsf] = y[dim + j][k] - y[dim + j][k - 1] -
                          dz * (y[j][k] + y[j][k - 1]) / 2;
      }
      // j = 2 * dim -> eta
      for (size_t n = 0; n < dim; n++)
      {
        s[dim + dim][0 * dim + indexv[n]]     = 0; // 9
        s[dim + dim][1 * dim + indexv[n]]     = 0; // 10
        s[dim + dim][2 * dim + indexv[n] + 1] = 0; // 14
        s[dim + dim][3 * dim + indexv[n] + 1] = 0; // 15
      }
      s[dim + dim][2 * dim]     = -1; // 13
      s[dim + dim][4 * dim + 1] = 1;  // 18
                                      //  Equations for E(k,k-1) for eta
      s[dim + dim][jsf] = y[dim + dim][k] - y[dim + dim][k - 1];
    }
  }
};