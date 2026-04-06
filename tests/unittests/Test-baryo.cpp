#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <memory>

using Approx = Catch::Approx;

#include <BSMPT/Kfactors/Kernels.h>
#include <BSMPT/baryo_fhck/BenchmarkModel.h>
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
  REQUIRE(transport.BAUeta.at(0).value() == Approx(-2.9673e-11).epsilon(1e-2));

  tmodel->VevProfile = Baryo::FHCK::VevProfileMode::FieldEquation;
  transport.Initialize();
  transport.SolveTransportEquation();
  REQUIRE(transport.BAUeta.at(0).value() == Approx(-3.87543e-11).epsilon(1e-2));
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

  transport.SolveTransportEquation();
  transport.CalculateBAU();
  CHECK(transport.bau == Approx(-3.87543e-11).epsilon(1e-2));
  for (auto &zi : tmodel->vacuumprofile->z)
    zi += 3 * tmodel->Lw;
  tmodel->vacuumprofile->GenerateSplines();

  transport.SolveTransportEquation();
  transport.CalculateBAU();
  CHECK(transport.bau == Approx(-3.87543e-11).epsilon(1e-2));
}

TEST_CASE("Test benchmark model for correctness.", "[BaryoBench]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  SetLogger({"--logginglevel::complete=true"});

  std::shared_ptr<BenchmarkModel> bmodel =
      std::make_shared<BenchmarkModel>(100., 0.1);
  
  TransportEquations transport(bmodel, 0.1);
  transport.SolveTransportEquation();
}