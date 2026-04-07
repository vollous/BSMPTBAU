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

TEST_CASE("Test baryo example_point_C2HDM", "[baryoFHCK]")
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
  Vacuum vac(
      0, 300, MinTracer, modelPointer, MultiStepPTMode::Default, 10, true);
  vac.setCoexPhases();

  std::shared_ptr<BSMPT::CoexPhases> coex =
      std::make_shared<BSMPT::CoexPhases>(vac.CoexPhasesList[0]);

  std::shared_ptr<TransportModel> tmodel =
      std::make_shared<TransportModel>(modelPointer,
                                       coex,
                                       0.01,
                                       coex->crit_temp,
                                       Baryo::FHCK::VevProfileMode::Kink);

  TransportEquations transport(tmodel, coex->crit_temp);

  tmodel->VevProfile = Baryo::FHCK::VevProfileMode::Kink;
  transport.SolveTransportEquation();
  REQUIRE(transport.BAUeta.at(0).value() == Approx(-2.93622e-11).epsilon(1e-2));

  tmodel->VevProfile = Baryo::FHCK::VevProfileMode::FieldEquation;
  transport.Initialize();
  transport.SolveTransportEquation();
  REQUIRE(transport.BAUeta.at(0).value() == Approx(-3.83627e-11).epsilon(1e-2));
}

TEST_CASE("Test baryo example_point_C2HDM z-invariance", "[baryoFHCK]")
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
  Vacuum vac(
      0, 300, MinTracer, modelPointer, MultiStepPTMode::Default, 10, true);
  vac.setCoexPhases();

  std::shared_ptr<BSMPT::CoexPhases> coex =
      std::make_shared<BSMPT::CoexPhases>(vac.CoexPhasesList[0]);

  std::shared_ptr<TransportModel> tmodel = std::make_shared<TransportModel>(
      modelPointer,
      coex,
      0.01,
      coex->crit_temp,
      Baryo::FHCK::VevProfileMode::FieldEquation);

  TransportEquations transport(tmodel, coex->crit_temp);

  transport.Initialize();
  transport.SolveTransportEquation();
  transport.CalculateBAU();
  CHECK(transport.bau == Approx(-3.83627e-11).epsilon(1e-2));
  for (auto &zi : tmodel->vacuumprofile->z)
    zi += 5 * tmodel->Lw;

  tmodel->vacuumprofile->GenerateSplines();

  // transport.Initialize();
  transport.SolveTransportEquation();
  transport.CalculateBAU();
  CHECK(transport.bau == Approx(-3.83627e-11).epsilon(1e-2));
}

TEST_CASE("Test baryo point. Difficult vacuum profile", "[baryoFHCK]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;

  const std::vector<double> example_point_C2HDM{
      /* lambda_1 = */ 1.2721750446204907,
      /* lambda_2 = */ 0.2483141342325688,
      /* lambda_3 = */ 8.574714969691337,
      /* lambda_4 = */ -3.209467738667644,
      /* Re(lambda_5) = */ -2.2755262749518543,
      /* Im(lambda_5) = */ -0.8030838099435237,
      /* Re(m_{12}^2) = */ 10780.786,
      /* tan(beta) = */ 17.147000,
      /* Yukawa Type = */ 1};

  using namespace BSMPT;
  SetLogger({"--logginglevel::complete=true"});
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  std::vector<double> TrueVacuum = {
      -3.256852516e-10, 14.33558551, 245.8304988, -0.001799560675};
  std::vector<double> FalseVacuum = {0, 0, 0, 0};
  double Tstar                    = 8.897363105;

  std::shared_ptr<TransportModel> tmodel = std::make_shared<TransportModel>(
      modelPointer, TrueVacuum, FalseVacuum, 0.1, Tstar);

  TransportEquations transport(tmodel, Tstar);

  transport.SolveTransportEquation();

  CHECK(transport.bau == Approx(-9.2609e-09).epsilon(1e-2));
}