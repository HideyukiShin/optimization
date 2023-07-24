#ifndef Q_newton_
#define Q_newton_

#include "../include/Eigen/Core"
#include "../include/Eigen/LU"
#include <utility>
#include <iostream>

template <size_t R>
class Q_newton 
{
  using X = Eigen::Matrix<double, R, 1>;
  using H = Eigen::Matrix<double, R, R>;

  X x;
  H hesse;
  double const epsilon = 1e-7;
  double const alpha = 0.001;

public:
  Q_newton(X const& init): x{init}, 
  hesse{Eigen::MatrixXd::Identity(R,R)} {}

  template<class F>
  X diff(F&& f, X const& xx)
  {
    X dx;
    X df;
    for(size_t i=0; i<xx.size(); ++i)
    {
      dx = xx;
      dx(i) += epsilon;
      df(i) = (f(dx) - f(xx))/epsilon;
    }
    return df;
  }

  template <class F>
  double linear_search(F&& f, X const& dir)
  {
    auto k=1;
    while(f(x+dir*k*alpha) <= f(x+dir*(k-1)*alpha))
      {
        ++k;
        if(k == 1000000)
        {
          break;
          throw std::runtime_error("linear search is failed");
        }
      }
      return (k-1)*alpha;
  }

  template <class F>
  std::pair<double, X> opt(F&& f)
  {
    auto steps = 0;
    X s;
    while(steps <= 100)
    {
      auto d = diff(f,x);
      if(d.norm() < 1e-5) return std::make_pair(f(x), x);
      auto delta_x = -hesse.inverse()*d;
     // auto a = linear_search(f, delta_x);
      auto next_x = x + delta_x;
      auto s = next_x - x;
      auto y = diff(f, next_x) - d;
      auto next_hesse = 
        hesse+(y*y.transpose())/(y.dot(s))
          -(hesse*s*s.transpose()*hesse)/(s.transpose()*hesse*s);
      x = next_x;
      hesse = next_hesse; 
      steps++;
    }
    return std::make_pair(f(x), x);
  }

};


#endif