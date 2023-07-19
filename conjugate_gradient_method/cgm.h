#ifndef CGM_
#define CGM_

#include "../include/Eigen/Core"
#include <utility>
#include <iostream>


template <size_t R>
class CGM
{
  using X = Eigen::Matrix<double, R, 1>;

  X x;
  double const epsilon = 0.00001;
  double const alpha = 0.001;

public:
  explicit CGM(X const& init) : x{init} {}

  template <class F>
  X diff(F&& f)
  {
    X dx = x;
    X df;
    for(size_t i=0; i<x.size(); ++i)
    {
      dx = x;
      dx(i) += epsilon;
      df(i) = (std::invoke(f, std::forward<X>(dx)) -
                  std::invoke(f, std::forward<X>(x)))/epsilon;
    }
    return df;
  }

  template <class F>
  double linear_search(F&& f, X const& s)
  {
    auto k = 1;
    while(std::invoke(f, std::forward<X>(x+k*alpha*s))
     <= std::invoke(f, std::forward<X>(x+(k-1)*alpha*s)))
     {
      ++k;
      if(k > 10000000000)
      {
        throw std::runtime_error("linear search is failed");
      }
     }
     return (k-1)*alpha;
  }

  template <class F>
  std::pair<double, X> opt(F&& f)
  {
    auto steps = 0;
    auto k=0;
    X s; 
    while(steps <= 100)
    {

      
      auto d = diff(std::forward<F>(f));
      if (d.norm() <= 0.1) 
      {
          return std::make_pair(
            std::invoke(f, std::forward<X>(x)), x);
      }
      if (k==0) s = -d;

      
      auto k = linear_search(std::forward<F>(f), s);
      x += k * s;
      if(k==R-1){ 
        k = 0;
        ++steps;
        continue;
      }

      auto d_next = -diff(std::forward<F>(f));
      auto beta = (d_next.dot(d_next))/d.dot(d);
      //auto beta = d_next.dot(d_next - d)/(d_next.dot(d_next));
      s = d_next + beta * s;
      ++k;
      ++steps;
    }

    return std::make_pair(
            std::invoke(f, std::forward<X>(x)), x);
  }

};


#endif
