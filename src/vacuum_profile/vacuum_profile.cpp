/**
 * @file
 */

#include <BSMPT/vacuum_profile/vacuum_profile.h>

namespace BSMPT
{
namespace VacuumProfileNS
{

VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t &dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In,
    // z knot list
    std::vector<double> &z_In,
    // field knot list
    std::vector<std::vector<double>> &path_In)
    : dim(dim_In)
    , TrueVacuum(TrueVacuum_In)
    , FalseVacuum(FalseVacuum_In)
    , V(V_In)
    , dV(dV_In)
    , Hessian(Hessian_In)
    , z(z_In)
    , path(path_In)
{
}

} // namespace VacuumProfileNS
} // namespace BSMPT