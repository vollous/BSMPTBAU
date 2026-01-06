#pragma once

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/SMparam.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/difeq_vacuum_profile.h>
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
} // namespace Baryo
} // namespace BSMPT
