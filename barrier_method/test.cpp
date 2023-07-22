
#include "../../include/Eigen/Core"
#include "Constraint.hpp"
#include"barrier_method.hpp"

double constraints_1(Eigen::Matrix<double, 2,1> a)
{
  return a.dot(a);
}

double constraint_2(Eigen::Matrix<double, 2, 1> a)
{
  return (a(0)-0.5)*((a(0)-0.5)) + (a(1)-0.7)*(a(1)-0.7);
}

double value_func(Eigen::Matrix<double, 2, 1> x)
{
  Eigen::Vector2d a;
  a << -1, -1;
  return   a.dot(x);
}

int main()
{
  auto init_pos =  Eigen::Vector2d{};
  init_pos << 0.1, 0;
  Constraint<2> c1(constraints_1, "<=", 1);
  Constraint<2> c2(constraint_2, ">=", 0.3);
  
  Barrier_method<2> b_method(init_pos);
  b_method.append_constraint(c1);
  b_method.append_constraint(c2);
 // std::cout <<   b_method.diff(value_func, init_pos) << "  " << b_method.cost_function(value_func, init_pos) << std::endl;
 auto result = b_method.opt(value_func);  std::cout << result.first << std::endl;  std::cout << result.second << std::endl;  std::cout << "penalty: " << constraints_1(result.second) << std::endl;
  std::cout << b_method.check_constraint(result.second, c1).first << std::endl;
  std::cout << b_method.check_constraint(result.second, c2).first << std::endl;
}