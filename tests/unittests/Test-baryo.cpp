#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/Kfactors/Kernels.h>
#include <BSMPT/baryo_fhck/TransportEquations.h>
#include <BSMPT/baryo_fhck/TransportModel.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/vacuum_profile.h>

TEST_CASE("Construct Kernel table", "[baryoKernels]")
{
  using namespace BSMPT::Baryo::FHCK;
  std::string path = "kernels/";
  double vw_step   = 0.064;

  for (int l = 21; l <= 30; l++)
  {
    std::cout << "Current moment: " << l << "\n\n";
    for (int type = 0; type <= 1; type++)
    {
      ParticleType PType =
          (type == 0 ? ParticleType::LeftFermion : ParticleType::Boson);
      std::string suffix = (type == 0 ? "_f.dat" : "_b.dat");
      std::cout << (type == 0 ? "Fermion" : "Boson") << "\n";
      Kernel Kern(l, 2);
      using std::chrono::duration;
      using std::chrono::duration_cast;
      using std::chrono::high_resolution_clock;
      using std::chrono::milliseconds;

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
      if (l != 1)
      {
        std::cout << "K-Kernel\n";
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
                   << Kern(KernelType::K, PType, x, vw) << "\n";
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
      if (l != 0)
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
        if (PType == ParticleType::LeftFermion)
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
      }
    }
  }
  exit(1);
  for (int type = 0; type <= 1; type++)
  {
    Kernel Kern(0, 0);
    ParticleType PType =
        (type == 0 ? ParticleType::LeftFermion : ParticleType::Boson);
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

TEST_CASE("Check example_point_C2HDM", "[baryoFHCK1]")
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

  // TODO remove this timing
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;

  auto t1 = high_resolution_clock::now();

  std::shared_ptr<TransportModel> tmodel =
      std::make_shared<TransportModel>(modelPointer,
                                       coex,
                                       0.01,
                                       coex->crit_temp,
                                       Baryo::FHCK::VevProfileMode::Kink);

  TransportEquations transport(tmodel, coex->crit_temp);

  transport.SolveTransportEquation();
  CHECK(transport.BAUEta.value() == Approx(-7.8868e-11).epsilon(1e-2));

  tmodel->VevProfile = Baryo::FHCK::VevProfileMode::FieldEquation;
  transport.Initialize();
  transport.SolveTransportEquation();
  CHECK(transport.BAUEta.value() == Approx(-1.18786e-10).epsilon(1e-2));

  auto t2 = high_resolution_clock::now();

  /* Getting number of milliseconds as an integer. */
  auto ms_int = duration_cast<milliseconds>(t2 - t1);

  /* Getting number of milliseconds as a double. */
  duration<double, std::milli> ms_double = t2 - t1;

  std::cout << ms_int.count() << "ms\n";
  std::cout << ms_double.count() << "ms\n";
}