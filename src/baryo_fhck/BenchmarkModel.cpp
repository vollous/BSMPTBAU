#include <BSMPT/baryo_fhck/BenchmarkModel.h>
#include <cmath>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

BenchmarkModel::BenchmarkModel(const double Tn_in,
                               const double vw_in,
                               const TruncationScheme &truncationscheme_in,
                               const double &truncationR_in)
    : TransportModel(nullptr,
                     {},
                     {},
                     vw_in,
                     Tn_in,
                     VevProfileMode::Kink,
                     truncationscheme_in,
                     truncationR_in)
    , vn(Tn_in)
    , wn(2. * Tn_in)
    , LAM(1000.)
{
  Lw = 5. / Tn_in;
}

BenchmarkModel::BenchmarkModel(const double vn_in,
                               const double wn_in,
                               const double Tn_in,
                               const double LAM_in,
                               const double Lw_in,
                               const double vw_in,
                               const TruncationScheme &truncationscheme_in,
                               const double &truncationR_in)
    : TransportModel(nullptr,
                     {},
                     {},
                     vw_in,
                     Tn_in,
                     VevProfileMode::Kink,
                     truncationscheme_in,
                     truncationR_in)
    , vn(vn_in)
    , wn(wn_in)
    , LAM(LAM_in)

{
  Lw = Lw_in;
  std::cout << "Instantiated\n";
}

void BenchmarkModel::Initialize()
{
  status = TransportModelStatus::Success;
  std::cout << "Initialized\n";
}

void BenchmarkModel::GenerateFermionMass(const std::vector<double> &zList,
                                         const bool &MakeTopMassPlot)
{
  (void)zList;
  (void)MakeTopMassPlot;
  std::cout << "FermionsMassesGenerated\n";
}

double BenchmarkModel::hvev(const double &z, const int &deriv)
{
  double res = 0.;
  if (deriv == 0)
    res += vn / 2. * (1 - tanh(z / Lw));
  else if (deriv == 1)
    res += -vn / (2. * Lw) / pow(cosh(z / Lw), 2);
  return res;
}

double BenchmarkModel::svev(const double &z, const int &deriv)
{
  double res = 0.;
  double u   = z / Lw;
  if (deriv == 0)
    res += wn / 2. * (1 + tanh(u));
  else if (deriv == 1)
    res += wn / (2. * Lw) / pow(cosh(u), 2);
  else if (deriv == 2)
    res += -wn / (Lw * Lw) * tanh(u) / pow(cosh(u), 2);
  return res;
}

void BenchmarkModel::GetFermionMass(const double &z,
                                    const size_t &fermion,
                                    double &m2,
                                    double &m2prime,
                                    double &thetaprime,
                                    double &theta2prime)
{
  // bottom quark has no mass in this model
  if (fermion == 2)
  {
    m2          = 0.;
    m2prime     = 0.;
    thetaprime  = 0.;
    theta2prime = 0.;
  }
  else
  {
    double h    = hvev(z, 0);
    double dh   = hvev(z, 1);
    double s    = svev(z, 0) / LAM;
    double ds   = svev(z, 1) / LAM;
    double d2s  = svev(z, 2) / LAM;
    double temp = sqrt(1. + s * s);
    m2          = yt * h * temp;
    m2prime     = yt * ((1. + s * s) * dh + h * s * ds) / temp;
    m2prime     = 2. * m2 * m2prime / (Tstar * Tstar);
    m2 *= m2 / (Tstar * Tstar);
    thetaprime  = ds / (1 + s * s);
    theta2prime = ((1 + s * s) * d2s - 2. * s * ds * ds) / pow(1. + s * s, 2);
  }
}

double BenchmarkModel::GetWMass(const double &z, const double &T)
{
  double h = hvev(z, 0);
  return h * g / 2. / T;
}

double BenchmarkModel::EWSBVEV(const double &z)
{
  return hvev(z, 0);
}

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT