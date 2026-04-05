#pragma once

#include <BSMPT/baryo_fhck/TransportModel.h>
#include <vector>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

class BenchmarkModel : public TransportModel
{
private:
  /**
   *@brief Input parameters of the benchmark model (see [2407.13639] Sec.7)
   *
   */
  double vn, wn, LAM;

  /**
   *@brief Top yukawa coupling and weak coupling constant
   */
  const double yt = 1., g = 0.65;

public:
  BenchmarkModel(
      const double Tn_in,
      const double vw_in,
      const TruncationScheme &truncationscheme_in = TruncationScheme::MinusVw,
      const double &truncationR_in                = 0);
  BenchmarkModel(
      const double vn_in,
      const double wn_in,
      const double Tn_in,
      const double LAM_in,
      const double Lw_in,
      const double vw_in,
      const TruncationScheme &truncationscheme_in = TruncationScheme::MinusVw,
      const double &truncationR_in                = 0);

  void Initialize() override;

  void GenerateFermionMass(const std::vector<double> &zList,
                           const bool &MakeTopMassPlot) override;

  double hvev(const double &z, const int &deriv);

  double svev(const double &z, const int &deriv);

  void GetFermionMass(const double &z,
                      const size_t &fermion,
                      double &m2,
                      double &m2prime,
                      double &thetaprime,
                      double &theta2prime) override;

  double GetWMass(const double &z, const double &T) override;

  double EWSBVEV(const double &z) override;
};
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT