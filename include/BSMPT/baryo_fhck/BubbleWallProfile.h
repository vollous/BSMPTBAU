#pragma once

#include <BSMPT/baryo_fhck/difeq_tunneling_path.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/SMparam.h>
#include <BSMPT/utility/utility.h>
#include <iostream>
#include <string>
#include <vector>

namespace BSMPT
{
namespace Baryo
{
void CalculateBubbleWallProfile(
    std::vector<std::vector<double>> &path,
    const std::function<std::vector<double>(std::vector<double>)> &dV);
}
} // namespace BSMPT
