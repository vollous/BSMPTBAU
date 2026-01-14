#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/Kfactors/Kernels.h>
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
  double vw_step   = 0.064;

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
        for (double i = 0.; i < 7.01; i += 0.1)
        {
          double x = -5 + i;
          x        = pow(10, x);
          if (l == 0)
            file << x << "\t" << Kern(KernelType::D, PType, x, 0.1) << "\n";
          else
            for (double j = 0.; j < 5.01; j += vw_step)
            {
              double vw = -5. + j;
              vw        = pow(10, vw);
              file << x << "\t" << vw << "\t"
                   << Kern(KernelType::D, PType, x, vw) << "\n";
            }
        }
        file.close();
      }
      std::cout << "K-Kernel\n";
      {
        std::string str = path + "K" + std::to_string(l) + suffix;
        std::ofstream file(str);
        for (double i = 0.; i < 7.01; i += 0.1)
        {
          double x = -5 + i;
          x        = pow(10, x);
          if (l == 0)
            file << x << "\t" << Kern(KernelType::K, PType, x, 0.1) << "\n";
          else if (l >= 2)
            for (double j = 0.; j < 5.01; j += vw_step)
            {
              double vw = -5. + j;
              vw        = pow(10, vw);
              file << x << "\t" << vw << "\t"
                   << Kern(KernelType::D, PType, x, vw) << "\n";
            }
        }
        file.close();
      }
      if (l == 0)
      {
        std::cout << "K4FH-Kernel\n";
        std::string str = path + "K4FH" + suffix;
        std::ofstream file(str);
        for (double i = 0.; i < 7.01; i += 0.1)
        {
          double x = -5 + i;
          x        = pow(10, x);
          if (l == 0)
            file << x << "\t" << Kern(KernelType::K4FH, PType, x, 0.1) << "\n";
        }
        file.close();
      }
      /* if (l != 0)
      {
        std::cout << "Q-Kernel\n";
        {
          std::string str = path + "Q" + std::to_string(l) + suffix;
          std::ofstream file(str);
          for (double i = 0.; i < 7.01; i += 0.1)
          {
            double x = -5 + i;
            x        = pow(10, x);
            for (double j = 0.; j < 5.01; j += vw_step)
            {
              double vw = -5. + j;
              vw        = pow(10, vw);
              file << x << "\t" << vw << "\t"
                   << Kern(KernelType::Q, PType, x, vw) << "\n";
            }
          }
          file.close();
        }
        if (PType == fermion)
        {
          std::cout << "Q8o-Kernel\n";
          {
            std::string str = path + "Q8o" + std::to_string(l) + suffix;
            std::ofstream file(str);
            for (double i = 0.; i < 7.01; i += 0.1)
            {
              double x = -5 + i;
              x        = pow(10, x);
              for (double j = 0.; j < 5.01; j += vw_step)
              {
                double vw = -5. + j;
                vw        = pow(10, vw);
                file << x << "\t" << vw << "\t"
                     << Kern(KernelType::Q8o, PType, x, vw) << "\n";
              }
            }
            file.close();
          }
          std::cout << "Q9o-Kernel\n";
          {
            std::string str = path + "Q9o" + std::to_string(l) + suffix;
            std::ofstream file(str);
            for (double i = 0.; i < 7.01; i += 0.1)
            {
              double x = -5 + i;
              x        = pow(10, x);
              for (double j = 0.; j < 5.01; j += vw_step)
              {
                double vw = -5. + j;
                vw        = pow(10, vw);
                file << x << "\t" << vw << "\t"
                     << Kern(KernelType::Q9o, PType, x, vw) << "\n";
              }
            }
            file.close();
          }
        }
      } */
    }
  }
  for (int type = 0; type <= 1; type++)
  {
    Kernel Kern(0, 0);
    ParticleType PType = (type == 0 ? Fermion : Boson);
    std::string suffix = (type == 0 ? "_f.dat" : "_b.dat");
    std::cout << (type == 0 ? "Fermion" : "Boson") << "\n";
    std::cout << "Rbar-Kernel\n";
    {
      std::string str = path + "Rbar" + suffix;
      std::ofstream file(str);
      for (double i = 0.; i < 7.01; i += 0.1)
      {
        double x = -5 + i;
        x        = pow(10, x);
        for (double j = 0.; j < 5.01; j += vw_step)
        {
          double vw = -5. + j;
          vw        = pow(10, vw);
          file << x << "\t" << vw << "\t" << Kern(KernelType::Rb, PType, x, vw)
               << "\n";
        }
      }
      file.close();
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


  TransportEquations transport(
      modelPointer,
      coex,
      0.01,
      coex->crit_temp,
      Baryo::FHCK::VevProfileMode::Kink); // TODO: rename this

  transport.SolveTransportEquation();

  REQUIRE(transport.BAUEta.value() == Approx(1.261614739e-10).epsilon(1e-4));

  transport.VevProfile = Baryo::FHCK::VevProfileMode::TunnelPath;
  transport.Initialize();
  transport.SolveTransportEquation();
  REQUIRE(transport.BAUEta.value() == Approx(9.895749505e-09).epsilon(1e-4));

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
  double Tstar                    = 143.7102777;
  double eps                      = 0.01;
  size_t NumberOfSteps            = 1000;
  std::vector<double> FalseVacuum = {0, 0, 0, 0};
  std::vector<double> TrueVacuum  = {0, 49.3578, 194.183, -1.30873};
  double itmax                    = 30;
  double conv                     = 1e-10;
  double slowc                    = 1e-2;
  int NB                          = 4;
  std::vector<double> zList(NumberOfSteps);
  MatDoub y(dim * 2, zList.size(), 0.);
  y.zero();
  VecDoub scalv(zList.size(), 246.22);
  VecInt indexv(2 * dim);
  VacuumProfileNS::ProfileSolverMode mode =
      VacuumProfileNS::ProfileSolverMode::Deriv;
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
    const double dfac = 10 * pow(cosh(10 * zList[i]), -2) / 2.;
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

TEST_CASE("Domain Wall lambda^4", "[baryoFHCKdomain]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

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

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);

  vacuumprofile.CalculateProfile();

  REQUIRE(1 == 1);
}
TEST_CASE("Bubble profile lambda^4", "[baryoFHCKdomain]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

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

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);
  vacuumprofile.CalculateProfile();

  REQUIRE(1 == 1);
}

TEST_CASE("Test indexv", "[baryoFHCK]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  size_t dim = 4;
  std::vector<double> TrueVacuum(dim, -1);
  std::vector<double> FalseVacuum(dim, 1);

  std::function<double(std::vector<double>)> V;
  std::function<std::vector<double>(std::vector<double>)> dV;
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian, 1);

  VecInt indexv, indexvField(8), indexvDeriv(8);

  indexvDeriv[0] = 0;
  indexvDeriv[1] = 1;
  indexvDeriv[2] = 2;
  indexvDeriv[3] = 3;
  indexvDeriv[4] = 4;
  indexvDeriv[5] = 5;
  indexvDeriv[6] = 6;
  indexvDeriv[7] = 7;

  indexvField[0] = 4;
  indexvField[1] = 5;
  indexvField[2] = 6;
  indexvField[3] = 7;
  indexvField[4] = 0;
  indexvField[5] = 1;
  indexvField[6] = 2;
  indexvField[7] = 3;

  vacuumprofile.mode = VacuumProfileNS::ProfileSolverMode::Field;
  indexv             = vacuumprofile.Calcindexv();
  for (size_t i = 0; i < 2 * dim; i++)
    CHECK(indexv[i] == indexvField[i]);

  vacuumprofile.mode = VacuumProfileNS::ProfileSolverMode::Deriv;
  indexv             = vacuumprofile.Calcindexv();
  for (size_t i = 0; i < 2 * dim; i++)
    CHECK(indexv[i] == indexvDeriv[i]);
}

TEST_CASE("Test center path", "[baryoFHCK]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

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

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);

  vacuumprofile.CalculateProfile();

  double center;
  vacuumprofile.CenterPath(center); // Center path
  vacuumprofile.CenterPath(center); // Recenter does nothing

  REQUIRE(center == Approx(0.).epsilon(1e-8));
}