#ifndef BARRIER_
#define BARRIER_

#include "../..//include/Eigen/Core"
#include "../../include/Eigen/LU"
#include "Constraint.hpp"
#include <utility>
#include <iostream>
#include <string>
#include <vector>

template <size_t R>
class Barrier_method
{
  enum class algorithm {SDM, CGM, NEWTON};
  using X = Eigen::Matrix<double, R, 1>;
  using Hesse = Eigen::Matrix<double, R, R>;

  X x;
  Hesse H;
  algorithm alg;

  std::vector<Constraint<R>> constraints;
  double const epsilon = 1e-7;
  double const alpha = 1e-5;
  double const eta = 1e-3;
  double weight = 1;
  double WEIGHT_LIMIT = 1000000;
  double STEP_LIMIT = 1000;
  double LINEAR_LIMIT = 1000000;

public:
  Barrier_method(X const& init, std::string const& a="sdm"):
    x{init}, H{}, alg(algorithm::SDM), constraints{}
    {
      if(a=="sdm" || a=="SDM") alg = algorithm::SDM;
      else if(a=="cgm"|| a=="CGM") alg = algorithm::CGM;
      else if(a=="newton" || a=="NEWTON") alg=algorithm::NEWTON;
    }

  void append_constraint(Constraint<R> const& c)
  {
    constraints.push_back(c);
  }

  
  std::pair<bool, double> check_constraint(X const& x_, 
  Constraint<R> c)
  {
    bool result_bool ; 
    double result_value = c.left(x_) - c.right;;
    switch(c.op)
    {
      case Constraint<R>::type::le: 
        result_bool = result_value <= 0;
        return std::make_pair(result_bool, -result_value);
      case Constraint<R>::type::lt:
        result_bool = result_value < 0;
        return std::make_pair(result_bool, -result_value);
      case Constraint<R>::type::ge: 
        result_bool = result_value > 0;
        return std::make_pair(result_bool, result_value);
      case Constraint<R>::type::gt:
        result_bool = result_value >= 0;
        return std::make_pair(result_bool, result_value);
      case Constraint<R>::type::eq:
        result_bool = result_value == 0;
        return std::make_pair(result_bool, result_value);
    }
    throw std::runtime_error("type is invalid!");
  }

  
  double barrier_function(double const& x_)
  {
    return 1.0/(x_*x_);
  }

  template <class F>
  double cost_function(F&& f, X const& x_)
  {
    double sum = 0;
    for(const auto& c: constraints)
    {
      if(c.op == Constraint<R>::type::eq)
      {
        throw std::runtime_error(
            "barrier method cannot apply equation constraint");
      }
      auto bv = check_constraint(x_, c);
      if(bv.first) sum += barrier_function(bv.second)/weight;
      else 
      {
        sum =  std::numeric_limits<double>::infinity();
        break;
      }
    }
    return f(x_) + sum;
  }

  template <class F>
  X diff(F&& f, X const& x_)
  {
    X df;
    for(auto i=0; i<x_.size(); ++i)
    {
      auto dx = x_;
      dx(i) += epsilon;
      df(i) = (cost_function(f, dx)-cost_function(f, x_))/epsilon;
    }
    return df;
  }

  template <class F>
  double linear_search(F&& f, X const& dir)
  {
    auto k = 1;
    while(cost_function(f, x+k*alpha*dir)
      <=cost_function(f, x+(k-1)*alpha*dir))
      {
        ++k;
        if(k>LINEAR_LIMIT) throw std::runtime_error("Linear search is failed");
      }
      return (k-1)*alpha;
  }

  template <class F>
  std::pair<double, X> opt(F&& f)
  {
    auto prev_cost = std::numeric_limits<double>::infinity();
    while(weight<WEIGHT_LIMIT)
    {
      auto steps=0;
      while(steps<=STEP_LIMIT)
      {
        auto d = -diff(f, x);
        if(d.norm() < 1e-9) {
          std::cout << "reached!" << std::endl;
          break;
        }
        auto a = linear_search(f, d);
        x += a*d;
        steps++;
        if(std::fabs(prev_cost-cost_function(f,x))<1e-9) break;
        prev_cost = cost_function(f,x);
        
      }
      weight += 1000;
    }
    return std::make_pair(f(x), x);
  }

};



#endif