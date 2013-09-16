// This file is part of PG.
//
// PG is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PG is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with PG.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

// include
// roboptim
#include <roboptim/core/differentiable-function.hh>

namespace pg
{

template<typename Type, int Size>
class AutoDiffFunction : public roboptim::DifferentiableFunction
{
public:
  typedef typename Type::scalar_t scalar_t;
  typedef Eigen::Matrix<scalar_t, Size, 1> result_ad_t;

public:
  AutoDiffFunction(int pbSize, int size, const std::string& name)
    : roboptim::DifferentiableFunction(pbSize, size, name)
    , res_(size)
  {
    res_.fill(scalar_t(0., Eigen::VectorXd::Zero(pbSize)));
  }
  ~AutoDiffFunction() throw()
  {}

  virtual void impl_compute(result_ad_t& res, const argument_t& x) const = 0;


  void impl_compute(result_t& res, const argument_t& x) const throw()
  {
    impl_compute(res_, x);

    for(int i = 0; i < outputSize(); ++i)
    {
      /// @todo Put in Type some function to extract value
      res[i] = res_[i].value();
    }
  }


  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
  {
    impl_compute(res_, x);

    for(int i = 0; i < outputSize(); ++i)
    {
      /// @todo Put in Type some function to extract derivatives
      jac.row(i) = res_[i].derivatives().transpose();
    }
  }


  void impl_gradient(gradient_t& gradient,
      const argument_t& x, size_type functionId) const throw()
  {
    impl_compute(res_, x);

    gradient = res_[functionId].derivatives().transpose();
  }

private:
  mutable result_ad_t res_;
};

} // namespace pg



