#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/Kfactors/Kernels.h>
#include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/baryo_fhck/TransportEquations.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/vacuum_profile.h>

TEST_CASE("Construct Kernel table", "[baryoKernels]")
{
  using namespace BSMPT;
  std::string path = "kernels/";

  for (int l = 0; l <= 2; l++)
  {
    std::cout << "Current moment: " << l << "\n\n";
    for (int type = 0; type <= 1; type++)
    {
      ParticleType PType = (type == 0 ? Fermion : Boson);
      std::string suffix = (type == 0 ? "_f.dat" : "_b.dat");
      std::cout << (type == 0 ? "Fermion" : "Boson") << "\n";
      Kernel Kern(l, 2);

      std::cout << "D-Kernel\n";
      {
        std::string str = path + "D" + std::to_string(l) + suffix;
        std::ofstream file(str);
        for (double i = 0.; i < 10.01; i += 0.1)
        {
          double x = -8 + i;
          x        = pow(10, x);
          if (l == 0)
            file << x << "\t" << Kern(KernelType::D, PType, x, 0.1) << "\n";
          else
            for (double j = 0.; j < 5.01; j += 0.064)
            {
              double vw = -5. + j;
              vw        = pow(10, vw);
              file << x << "\t" << vw << "\t"
                   << Kern(KernelType::D, PType, x, vw) << "\n";
            }
        };
        file.close();
      }
      if (l != 0)
      {
        /* std::cout << "Q-Kernel\n";
        {
          std::string str = path + "Q" + std::to_string(l) + suffix;
          std::ofstream file(str);
          for (double i = 0.; i < 10.01; i += 0.1)
          {
            double x = -8 + i;
            x        = pow(10, x);
            file << x << "\t" << Kern(KernelType::Q, PType, x, vw) << "\n";
          };
          file.close();
        }
        std::cout << "Q8o-Kernel\n";
        {
          std::string str = path + "Q8o" + std::to_string(l) + suffix;
          std::ofstream file(str);
          for (double i = 0.; i < 10.01; i += 0.1)
          {
            double x = -8 + i;
            x        = pow(10, x);
            file << x << "\t" << Kern(KernelType::Q8o, PType, x, vw) << "\n";
          };
          file.close();
        } */
        /* std::cout << "Q9o-Kernel\n";
        {
          std::string str = path + "Q9o" + std::to_string(l) + suffix;
          std::ofstream file(str);
          for (double i = 0.; i < 10.01; i += 0.1)
          {
            double x = -8 + i;
            x        = pow(10, x);
            std::cout << x << "\n";
            file << x << "\t" << Kern(KernelType::Q9o, PType, x, vw) << "\n";
          };
          file.close();
        } */
      }
    }
  }

  exit(0);

  REQUIRE(1 == 1);
}

TEST_CASE("Check example_point_C2HDM", "[baryoFHCK]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;

  const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                                /* lambda_2 = */ 0.274365,
                                                /* lambda_3 = */ 4.71019,
                                                /* lambda_4 = */ -2.23056,
                                                /* Re(lambda_5) = */ -2.43487,
                                                /* Im(lambda_5) = */ 0.124948,
                                                /* Re(m_{12}^2) = */ 2706.86,
                                                /* tan(beta) = */ 4.64487,
                                                /* Yukawa Type = */ 1};

  using namespace BSMPT;
  SetLogger({"--logginglevel::complete=true"});
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  Vacuum vac(0, 300, MinTracer, modelPointer, 1, 10, true);
  vac.setCoexPhases();

  std::shared_ptr<BSMPT::CoexPhases> coex =
      std::make_shared<BSMPT::CoexPhases>(vac.CoexPhasesList[0]);

  TransportEquations transport(modelPointer, coex, 0.01, coex->crit_temp);

  transport.SolveTransportEquation();

  exit(0);

  REQUIRE(1 == 1);
}

TEST_CASE("VEV Profile", "[baryo123]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;

  const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                                /* lambda_2 = */ 0.274365,
                                                /* lambda_3 = */ 4.71019,
                                                /* lambda_4 = */ -2.23056,
                                                /* Re(lambda_5) = */ -2.43487,
                                                /* Im(lambda_5) = */ 0.124948,
                                                /* Re(m_{12}^2) = */ 2706.86,
                                                /* tan(beta) = */ 4.64487,
                                                /* Yukawa Type = */ 1};

  using namespace BSMPT;
  SetLogger({"--logginglevel::complete=true"});
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  size_t dim                      = 4;
  double Tstar                    = 143.71;
  double eps                      = 0.01;
  size_t NumberOfSteps            = 1000;
  std::vector<double> FalseVacuum = {0, 0, 0, 0};
  std::vector<double> TrueVacuum  = {0, 49.3578, 194.183, -1.30873};
  double itmax                    = 30;
  double conv                     = 1e-10;
  double slowc                    = 1e-2;
  int NB                          = 4;
  VecDoub zList(NumberOfSteps);
  MatDoub y(dim * 2, zList.size(), 0.);
  y.zero();
  VecDoub scalv(zList.size(), 246.22);
  VecInt indexv(2 * dim);
  VacuumProfileNS::ProfileSolverMode mode =
      VacuumProfileNS::ProfileSolverMode::Field;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
  {
    // Potential wrapper
    std::vector<double> res = modelPointer->MinimizeOrderVEV(vev);
    return modelPointer->VEff(res, Tstar);
  };

  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg) { return NablaNumerical(arg, V, eps); };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg) { return HessianNumerical(arg, V, eps); };

  for (size_t i = 0; i < NumberOfSteps; i++)
  {
    double temp = ((i - (NumberOfSteps - 1) / 2.)) / ((NumberOfSteps - 1) / 2.);
    zList[i]    = 1 * temp;

    const double fac  = (tanh(10 * zList[i]) + 1.) / 2.;
    const double dfac = pow(cosh(zList[i]), -2) / 2.;
    for (size_t d = 0; d < dim; d++)
    {
      y[d][i]       = FalseVacuum[d] * dfac + TrueVacuum[d] * (-dfac);
      y[dim + d][i] = FalseVacuum[d] * fac + TrueVacuum[d] * (1 - fac);
    }
  }

  VacuumProfileNS::Difeq_VacuumProfile difeq_vacuumprofile(
      mode, dim, zList, TrueVacuum, FalseVacuum, V, dV, Hessian);

  for (int i = 0; i < 2 * dim; i++)
    indexv[i] = i;

  RelaxOde solvde(
      itmax, conv, slowc, scalv, indexv, NB, y, difeq_vacuumprofile);

  exit(0);

  REQUIRE(1 == 1);
}

TEST_CASE("Domain Wall lambda^4 Mode=Deriv", "[baryoFHCKdomain]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  double eps                      = 0.01;
  size_t NumberOfSteps            = 10000;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};
  double itmax                    = 20;
  double conv                     = 1e-8;
  double slowc                    = 1;
  int NB                          = dim;
  ProfileSolverMode mode          = ProfileSolverMode::Deriv;
  VecDoub zList(NumberOfSteps);
  VecDoub scalv(2 * dim, 1);
  std::vector<double> z(NumberOfSteps);
  std::vector<std::vector<double>> path;
  VecInt indexv(2 * dim + 1);
  for (size_t i = 0; i < 2 * dim + 1; i++)
  {
    // 0 1 2 3 4 ...
    indexv[i] = i;
    // switch index of field and deriv
    // 0 1 2 3 4 -> 2 3 0 1 4
    indexv[i] += dim * ((i < dim) * 2 - 1) *
                 /* dont update last element*/
                 (i < 2 * dim) *
                 /* reordeing only necessary for dirichlet */
                 (mode == ProfileSolverMode::Field);
  }

  double m   = 0;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  // Initial solution
  MatDoub y(dim * 2, zList.size(), 0.);
  y.zero();
  for (size_t i = 0; i < NumberOfSteps; i++)
  {
    double temp = ((i - (NumberOfSteps - 1) / 2.)) / ((NumberOfSteps - 1) / 2.);
    zList[i]    = 10 * temp;

    y[0][i] = sqrt(2) * pow(cosh(sqrt(2) * zList[i]), -2);
    y[1][i] = tanh(sqrt(2) * zList[i]) * (1 + sin(sqrt(2) * zList[i]) / 10.);
  }

  /* VacuumProfile::Difeq_VacuumProfile difeq_vacuumprofile(
      mode, dim, zList, TrueVacuum, FalseVacuum, V, dV, Hessian);
  RelaxOde solvde(
      itmax, conv, slowc, scalv, indexv, NB, y, difeq_vacuumprofile); */

  REQUIRE(1 == 1);
}

TEST_CASE("Domain Wall lambda^4 Mode=Field", "[baryoFHCKdomain]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  double eps                      = 0.01;
  size_t NumberOfSteps            = 10000;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};
  double itmax                    = 20;
  double conv                     = 1e-8;
  double slowc                    = 1;
  int NB                          = dim;
  VacuumProfileNS::ProfileSolverMode mode =
      VacuumProfileNS::ProfileSolverMode::Field;
  VecDoub zList(NumberOfSteps);
  VecDoub scalv(2 * dim, 1);

  VecInt indexv(2 * dim + 1);
  for (size_t i = 0; i < 2 * dim + 1; i++)
  {
    // 0 1 2 3 4 ...
    indexv[i] = i;
    // switch index of field and deriv
    // 0 1 2 3 4 -> 2 3 0 1 4
    indexv[i] += dim * ((i < dim) * 2 - 1) *
                 /* dont update last element*/
                 (i < 2 * dim) *
                 /* reordeing only necessary for dirichlet */
                 (mode == VacuumProfileNS::ProfileSolverMode::Field);
  }

  double m   = 0.1;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  // Initial solution
  MatDoub y(dim * 2, zList.size(), 0.);
  y.zero();
  for (size_t i = 0; i < NumberOfSteps; i++)
  {
    double temp = ((i - (NumberOfSteps - 1) / 2.)) / ((NumberOfSteps - 1) / 2.);
    zList[i]    = 10 * temp;

    y[0][i] = sqrt(2) * pow(cosh(sqrt(2) * zList[i]), -2);
    y[1][i] = tanh(sqrt(2) * zList[i]) * (1 + sin(sqrt(2) * zList[i]) / 10.);
  }

  VacuumProfileNS::Difeq_VacuumProfile difeq_vacuumprofile(
      mode, dim, zList, TrueVacuum, FalseVacuum, V, dV, Hessian);
  RelaxOde solvde(
      itmax, conv, slowc, scalv, indexv, NB, y, difeq_vacuumprofile);

  REQUIRE(1 == 1);
}
