#include <BSMPT/utility/NumericalIntegration.h>
#include <iomanip>

void save(const char *filename, const double &x, const std::vector<double> &y)
{
  std::cout << std::setprecision(16);
  std::ofstream outfile(filename, std::ios::out | std::ios::app);
  outfile.seekp(0, std::ios::end);
  outfile << x << "\t";
  for (auto it : y)
  {
    outfile << it << "\t";
  }
  outfile << "\n";
  outfile.close();
}

double setStep(double hnow, double err)
{
  double hnext;
  if (err > 1)
  {
    hnext = 0.9 * 1 / (sqrt(sqrt(err))) * hnow;
    if (hnext < 0.1 * hnow)
    {
      hnext = 0.1 * hnow;
    }
  }
  else
  {
    if (err > 0.0006)
    {
      hnext = 0.9 * std::pow(err, -0.2) * hnow;
    }
    else
    {
      hnext = 5 * hnow;
    }
  }
  return hnext;
}