/**
 * @file
 */

#include <BSMPT/baryo_fhck/BubbleWallProfile.h>

namespace BSMPT
{
namespace Baryo
{
void CalculateBubbleWallProfile(
    std::vector<std::vector<double>> &path,
    const std::function<std::vector<double>(std::vector<double>)> &dV)
{
  std::cout << "Working\n" << path << "\n\n\n";
  for (auto point : path)
    std::cout << "dV(" << point << ")"
              << " = " << dV(point) << "\n";
  exit(0);
}

} // namespace Baryo
} // namespace BSMPT