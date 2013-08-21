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
  StdCostFunc(PGData<Type>* pgdata, std::vector<std::vector<double>> q)
    : parent_t(pgdata->pbSize(), 1, "StdCostFunc")
    , pgdata_(pgdata)
    , tq_(std::move(q))
  {}


  void impl_compute(result_ad_t& res, const argument_t& x) const throw()
  {
    pgdata_->x(x);

    // compute posture task
    scalar_t posture = 0.;
    const std::vector<std::vector<scalar_t>>& q = pgdata_->q();
    for(int i = 0; i < pgdata_->multibody().nrJoints(); ++i)
    {
      if(pgdata_->multibody().joint(i).params() == 1)
      {
        posture += std::pow(tq_[i][0] - q[i][0], 2);
      }
    }

    res(0) = posture;
  }

private:
  PGData<Type>* pgdata_;
  std::vector<std::vector<double>> tq_;
};

} // namespace pg
