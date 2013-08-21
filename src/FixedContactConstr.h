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

// PG
#include "PGData.h"

namespace pg
{

template<typename Type>
class FixedContactConstr : public roboptim::DifferentiableFunction
{
public:
  typedef typename Type::scalar_t scalar_t;

public:
  FixedContactConstr(PGData<Type>* pgdata, int bodyId, const Eigen::Vector3d& pos)
    : roboptim::DifferentiableFunction(pgdata->pbSize(), 3, "FixedContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , pos_(pos.cast<scalar_t>())
  {}
  ~FixedContactConstr() throw()
  { }

  void impl_compute(result_t& res, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    Eigen::Vector3<scalar_t> resDiff = pgdata_->fk().bodyPosW()[bodyIndex_].translation() - pos_;
    for(int i = 0; i < 3; ++i)
    {
      res[i] = resDiff[i].value();
    }
  }

  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    Eigen::Vector3<scalar_t> resDiff = pgdata_->fk().bodyPosW()[bodyIndex_].translation() - pos_;
    for(int i = 0; i < 3; ++i)
    {
      jac.row(i) = resDiff[i].derivatives().transpose();
    }
  }

  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* argument */, size_type /* functionId */) const throw()
  { }

private:
  PGData<Type>* pgdata_;

  int bodyIndex_;
  Eigen::Vector3<scalar_t> pos_;
};

} // namespace pg
