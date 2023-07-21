#include"Penalty_method.hpp"
#include "../../include/Eigen/Core"

double constraints_1(Eigen::Matrix<double, 2,1> a)
{
  return a.dot(a);
}

double value_func(Eigen::Matrix<double, 2, 1> x)
{
  Eigen::Vector2d a;
  a << -1, -1;
  return  a.dot(x);
}

int main()
{
  auto init_pos =  Eigen::Vector2d{};
  init_pos << 10, 0.5;
  Constraint<2> c1(constraints_1, "<=", 1);
  Penalty_method<2> p_method(init_pos);
  p_method.append_constraint(c1);

  auto result = p_method.opt(value_func);
  std::cout << result.first << std::endl;
  std::cout << result.second << std::endl;
  std::cout << "penalty: " << constraints_1(result.second) << std::endl;
}