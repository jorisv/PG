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
// PG
#include "AutoDiffFunction.h"
#include "PGData.h"

namespace pg
{

template<typename Type>
class StdCostFunc : public AutoDiffFunction<Type, 1>
{
public:
  typedef AutoDiffFunction<Type, 1> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  StdCostFunc(PGData<Type>* pgdata, std::vector<std::vector<double>> q,
              double postureScale, double torqueScale)
    : parent_t(pgdata->pbSize(), 1, "StdCostFunc")
    , pgdata_(pgdata)
    , tq_(std::move(q))
    , postureScale_(postureScale)
    , torqueScale_(torqueScale)
  {}


  void impl_compute(result_ad_t& res, const argument_t& x) const throw()
  {
    pgdata_->x(x);

    // compute posture task
    scalar_t posture = scalar_t(0., Eigen::VectorXd::Zero(this->inputSize()));
    if(postureScale_ > 0.)
    {
      const std::vector<std::vector<scalar_t>>& q = pgdata_->q();
      for(int i = 0; i < pgdata_->multibody().nrJoints(); ++i)
      {
        if(pgdata_->multibody().joint(i).params() == 1)
        {
          posture += std::pow(tq_[i][0] - q[i][0], 2);
        }
      }
    }

    // compute torque task
    scalar_t torque = scalar_t(0., Eigen::VectorXd::Zero(this->inputSize()));
    if(torqueScale_ > 0.)
    {
      const auto& torqueVec = pgdata_->id().torque();
      for(int i = 0; i < pgdata_->multibody().nrJoints(); ++i)
      {
        for(int j = 0; j < pgdata_->multibody().joint(i).dof(); ++j)
        {
          torque += std::pow(torqueVec[i](j), 2);
        }
      }
    }

    res(0) = posture*postureScale_ + torque*torqueScale_;
  }

private:
  PGData<Type>* pgdata_;
  std::vector<std::vector<double>> tq_;
  double postureScale_;
  double torqueScale_;
};

} // namespace pg
