#ifndef NEWTONS_
#define NEWTONS_

#include "../include/Eigen/Core"
#include "../include/Eigen/LU"
#include <utility>
#include <iostream>

template <size_t R>
class Newton_method
{
  using X = Eigen::Matrix<double, R, 1>;
  using H = Eigen::Matrix<double, R, R>;

  X x;
  double const epsilon = 0.0001;

public:
  Newton_method(X const& init): x{init} {}

  template <class F>
  double call_func(F&& f, X y)
  {
    if(std::is_class<F>{})
    {
      return f(y);
    }
    else 
    {
      return std::invoke(f, std::forward<X>(y));
    }
  }

  template <class F>
  X diff(F&& f, X const& xx)
  {
    X dx = xx;
    X df;
    for(size_t i=0; i<xx.size(); ++i)
    {
      dx = xx;
      dx(i) += epsilon;
      df(i) = (call_func(std::forward<F>(f), dx) 
        - call_func(std::forward<F>(f), xx))/epsilon; 
    }
    return df;
  }

  template <class F>
  H Hessian(F&& f)
  {
    H m;
    const auto df = diff(f, x);
    for(size_t i=0; i<R; ++i)
    {
      for(size_t j=0; j<R; ++j)
      {
        X dx = x;
        dx(j) += epsilon;
        auto df_ = diff(f,dx);
        m(i,j) = (df_(i)-df(i))/epsilon;
      }
    }
    return m;
  }

  template<class F>
  std::pair<double, X> opt(F&& f)
  {
    auto steps = 0;
    auto k = 0;
    X s;
    while(steps <= 100)
    {
      auto d = diff(f, x);
      if(d.norm() < 0.00001) 
        return std::make_pair(call_func(f, x), x);
      auto H = Hessian(f);
      auto delta_x = -H.inverse()*d;
      x += delta_x;
    }

    return std::make_pair(call_func(f, x), x);
  }

};

#endif