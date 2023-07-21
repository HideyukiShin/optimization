#ifndef PENALTY_
#define PENALTY_

#include "../../include/Eigen/Core"
#include "../../include/Eigen/LU"
#include <utility>
#include <iostream>
#include <string>
#include <vector>


template <size_t R>
class Penalty_method;

template <size_t R>
class Constraint
{
  
  using X = Eigen::Matrix<double, R, 1>;
  enum class type { eq, le, lt, ge, gt};

  std::function<double(X)> left;
  type op;
  double right;

public:
  Constraint(std::function<double(X)> const& left_, 
   std::string const& op_=nullptr, double const right_=0): 
    left{left_}, op{type::le}, right{right_}
    {
      if(op_ == "=")        op=type::eq;
      else if(op_ == "<")   op=type::lt;
      else if(op_ == "<=")  op=type::le;
      else if(op_ == ">")   op=type::gt;
      else if(op_ == ">=")  op=type::ge;
    }
    friend class Penalty_method<R>;
};

template <size_t R>
class Penalty_method
{
  enum class algorithm {SDM, CGM, NEWTON};
  using X = Eigen::Matrix<double, R, 1>;
  using Hesse = Eigen::Matrix<double, R, R>;
  
  X x;
  Hesse H;
  algorithm alg;

  //std::vector<std::vector<int>> constraints;
  std::vector<Constraint<R>> constraints;
  double const epsilon = 1e-7;
  const double alpha = 1e-5;
  const double eta = 1e-2;
  double weight = 1; 
  double WEIGHT_LIMIT = 100000;
  double STEP_LIMIT = 100000;

public:
  Penalty_method(X const& init, std::string const& a="sdm"):
    x{init}, H{}, alg{algorithm::SDM}, constraints{}
    {
      if(a == "sdm" || a=="SDM") alg = algorithm::SDM;
      else if(a == "cgm" || a=="CGM") alg = algorithm::CGM;
      else if(a =="newton" || a == "NEWTON") alg = algorithm::NEWTON;
    }

  void append_constraint(Constraint<R> const& c)
  {
    constraints.push_back(c);
  }

  template <class F>
  double penalty_function(F&& f, X const& x_)
  {
    double sum = 0;
    for(const auto& c: constraints)
    {
      if(c.op == Constraint<R>::type::eq){
        auto error = -c.left(x_) + c.right;
        auto k = error * error;
        sum += k;
      }
      else if(c.op == Constraint<R>::type::ge 
          || c.op == Constraint<R>::type::gt){
        auto error = -c.left(x_) + c.right;
        auto k = error >=0 ? error*error : 0;
        sum += k;
      }
      else if(c.op == Constraint<R>::type::le 
         || c.op == Constraint<R>::type::lt){
        auto error = c.left(x_) - c.right;
        auto k = error >=0 ? error*error : 0;
        sum += k;
        
      }
    }
    return f(x_) + weight * sum;
  }

  template <class F>
  X diff(F&& f, X const& xx)
  {
    X df;
    for(auto i=0; i<xx.size(); ++i)
    {
      auto dx = xx;
      dx(i) += epsilon;
      df(i) = (penalty_function(f,dx)-penalty_function(f,xx))/epsilon;
    }
    return df;
  }

    template <class F>
  double linear_search(F&& f, X const& s)
  {
    auto k = 1;
    while(penalty_function(f, x+k*alpha*s)
     <= penalty_function(f, x+(k-1)*alpha*s))
     {
      ++k;
      if(k > 10000000000)
      {
        throw std::runtime_error("linear search is failed");
      }
     }
     return (k-1)*alpha;
  }

  // TODO Implement other algorithm
  template <class F>
  std::pair<double, X> opt(F&& f)
  {
    auto prev_pena = std::numeric_limits<double>::infinity();
    while(weight<WEIGHT_LIMIT){
      auto steps=0;
      while(steps <= STEP_LIMIT)
      {
        auto d = -diff(f, x);
        if(d.norm() < 1e-5 ) {
          std::cout << "reached" << std::endl;
          break;
        }
        auto a = linear_search(f, d);
        x += a * d;
        steps++;
        if(std::fabs(prev_pena-penalty_function(f, x))<1e-9) break;
        prev_pena = penalty_function(f,x);
      }
    weight += 10;
    
  }
  return std::make_pair(f(x), x);
  }
}; 




#endif