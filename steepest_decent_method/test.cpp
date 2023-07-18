#include "../include/Eigen/Core"
#include <iostream>
#include "sdm.h"

double value_func(Eigen::Vector2d const& x)
{
  auto a = Eigen::Vector2d {};
  a(0) = 1;
  a(1) = 12;
  return x.dot(x) + a.dot(x) + 32;
}

template <class F, class R>
double value_(F&& f, R& init)
{
  return std::invoke(f, std::forward<R>(init));
}

int main(){
  
   auto init_pose = Eigen::Matrix<double, 2,1> {};
   SDM<2> sdm(init_pose, 0.001);
  // std::cout << value_func(init_pose) << std::endl;
   //std::cout << value_(value_func, init_pose) << std::endl;
  // std::cout << sdm.diff(value_func) << std::endl;
  //SDM<2> sdm(0.0, 0.001);


  auto result = sdm.opt(value_func);
  std::cout << result.first  << std::endl;
  std::cout << result.second << std::endl;

}

