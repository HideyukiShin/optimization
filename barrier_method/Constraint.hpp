#ifndef CONSTRAINTS_
#define CONSTRAINTS_

#include "../../include/Eigen/Core"
#include "../../include/Eigen/LU"
#include <utility>
#include <iostream>
#include <string>

template <size_t R>
struct Constraint
{
  using X = Eigen::Matrix<double, R, 1>;
  enum class type {eq, le, lt, ge, gt};
  std::function<double(X)> left;
  type op;
  double right;

  Constraint(std::function<double(X)> const& left_,
    std::string const& op_=nullptr, double const right_ = 0):
      left{left_}, op{type::le}, right{right_}
      {
        if(op_ == "=")        op=type::eq;
        else if(op_ == "<")   op=type::lt;
        else if(op_ == "<=")  op=type::le;
        else if(op_ == ">")   op=type::gt;
        else if(op_ == ">=")  op=type::ge;        
      }
};

#endif