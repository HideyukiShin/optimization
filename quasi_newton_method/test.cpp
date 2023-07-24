#include "../include/Eigen/Core"
#include <iostream>
#include "quasi_newton.hpp"

double value_func(Eigen::Vector3d const& x)
{
  auto m = Eigen::Matrix<double, 3,3> {};
  m << 1, 4,4,5,4,6,7,9,10;

  auto a = Eigen::Vector3d {};
  a(0) = 1;
  a(1) = 12;
  a(2)=2;
  
  return x.transpose()*m*x + a.dot(x) + 32;
}

Eigen::Vector2d solve_equation(Eigen::Matrix2d const& A, 
  Eigen::Vector2d const& b)
  {
    auto func = [&A, &b](Eigen::Vector2d const& x)
      {return static_cast<double>(x.transpose()*A*x)/2 - b.dot(x) ;};
    auto y = Eigen::Vector2d{};
    y << 1,1;
    Q_newton<2> q(y);
   // std::cout << q.diff(func, y);
    auto result =  q.opt(func);
     std::cout << "result  " <<  result.second << std::endl; 
     return result.second;
    //return y;
  }

int main()
{
  auto init_pose = Eigen::VectorXd::Random(3);
   Q_newton<3> newton(init_pose);

  auto result = newton.opt(value_func);

  std::cout << result.first << std::endl;
  std::cout << result.second << std::endl;

  Eigen::Matrix2d m;
  Eigen::Vector2d b;
  m << 1,3,3,1;
  b << 1,0;
  auto s = solve_equation(m, b);
  std::cout << m * s << std::endl;
}