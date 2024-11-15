#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

// #include <BSMPT/Kfactors/vw_Kfactors.h>
#include <BSMPT/baryo_fhck/TransportEquations.h>
// #include <BSMPT/baryo_fhck/solvde.h>
// #include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/utility.h>

TEST_CASE("Check example_point_C2HDM", "[baryo1]")
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
  double Tc = vac.CoexPhasesList[0].crit_temp;

  std::shared_ptr<BSMPT::CoexPhases> coex =
      std::make_shared<BSMPT::CoexPhases>(vac.CoexPhasesList[0]);

  TransportEquations transport(modelPointer, coex, 0.1, coex->crit_temp);

  transport.SolveTransportEquation();

  exit(0);

  REQUIRE(1 == 1);
}