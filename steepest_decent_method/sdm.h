#ifndef SDM_
#define SDM_

#include "../include/Eigen/Core"
#include <utility>
#include <cstdlib>

template <size_t R>
class SDM
{
  using X = Eigen::Matrix<double, R, 1>;
  X  x;
  const double eta;
  const double epsilon = 0.000001;

public:
  SDM(X const& x_, double const& eta_);
  
  template <class F>
  X diff(F&& f);

  template <class F>
  std::pair<double, X> opt(F&& f);
};


#endif 