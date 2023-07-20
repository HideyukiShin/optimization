#include "../include/Eigen/Core"
#include <iostream>
#include "newton.hpp"

double value_func(Eigen::Vector4d const& x)
{
  auto a = Eigen::Vector4d {};
  a(0) = 1;
  a(1) = 12;
  a(2) = -4;
  a(3) = 6;
  return x.dot(x) + a.dot(x) + 32;
}

int main()
{
  auto init_pose = Eigen::Matrix<double, 4,1> {};
  Newton_method<4> newton(init_pose);

  auto result = newton.opt(value_func);

  std::cout << result.first << std::endl;
  std::cout << result.second << std::endl;
}