#include "sdm.h"
#include "../include/Eigen/Core"
#include <utility>
#include <iostream>

template <size_t R>
 SDM<R>::SDM(SDM::X const& x_, double const& eta_):
  x{x_}, eta{eta_} {}


  template <size_t R>
  template <class F>
  std::pair<double, typename SDM<R>::X>  SDM<R>::opt(F&& f)
  {
    auto steps = 0;
    while(steps<=10000000)
    {
      auto df = diff(std::forward<F>(f));
      if (df.norm() <= 0.000001) 
          return std::make_pair(std::invoke(f, std::forward<X>(x)), x);
      x -= eta*df;
      ++steps;
     // std::cout << x << std::endl;
    }
    return std::make_pair(std::invoke(f, std::forward<X>(x)), x);
  }




  template <size_t R>
  template <class F>
  SDM<R>::X SDM<R>::diff(F&& f)
  {
      X dx = x;
      X df ;
      for(auto i=0; i < x.size(); ++i)
      {
         dx = x;
        dx(i) += epsilon;
         df(i) = (std::invoke(f, std::forward<X>(dx)) - std::invoke(f, 
         std::forward<X>(x)))/epsilon;
      }
      
    //  return X{1,2};
      return df;
  }
