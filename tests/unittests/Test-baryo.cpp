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

  tmodel->VevProfile = Baryo::FHCK::VevProfileMode::Kink;
  transport.SolveTransportEquation();
  CHECK(transport.BAUeta.at(0).value() == Approx(-9.08465e-11).epsilon(1e-2));

  tmodel->VevProfile = Baryo::FHCK::VevProfileMode::FieldEquation;
  transport.Initialize();
  transport.SolveTransportEquation();
  CHECK(transport.BAUeta.at(0).value() == Approx(-1.36147e-10).epsilon(1e-2));

  auto t2 = high_resolution_clock::now();

  /* Getting number of milliseconds as an integer. */
  auto ms_int = duration_cast<milliseconds>(t2 - t1);

  /* Getting number of milliseconds as a double. */
  duration<double, std::milli> ms_double = t2 - t1;

  std::cout << ms_int.count() << "ms\n";
  std::cout << ms_double.count() << "ms\n";
}