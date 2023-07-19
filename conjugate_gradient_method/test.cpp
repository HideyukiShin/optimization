#include "../include/Eigen/Core"
#include <iostream>
#include "cgm.h"

double value_func(Eigen::Vector4d const& x)
{
  auto a = Eigen::Vector4d {};
  a(0) = 1;
  a(1) = 12;
  a(2) = -4;
  a(3) = 6;
  return x.dot(x) + a.dot(x) + 32;
}


int main(){
  
   auto init_pose = Eigen::Matrix<double, 4,1> {};
   CGM<4> cgm(init_pose);
  // std::cout << value_func(init_pose) << std::endl;
   //std::cout << value_(value_func, init_pose) << std::endl;
  // std::cout << sdm.diff(value_func) << std::endl;
  //SDM<2> sdm(0.0, 0.001);


  auto result = cgm.opt(value_func);
  std::cout << result.first  << std::endl;
  std::cout << result.second << std::endl;

}
