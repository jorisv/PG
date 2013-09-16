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
class FrictionConeConstr : public AutoDiffFunction<Type, Eigen::Dynamic>
{
public:
  typedef AutoDiffFunction<Type, Eigen::Dynamic> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  FrictionConeConstr(PGData<Type>* pgdata)
    : parent_t(pgdata, pgdata->pbSize(), pgdata->nrForcePoints(), "FrictionConeConstr")
    , pgdata_(pgdata)
  {}
  ~FrictionConeConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
    int i = 0;
    for(const auto& fd: pgdata_->forceDatas())
    {
      double muS = std::pow(fd.mu, 2);
      for(std::size_t fi = 0; fi < fd.points.size(); ++fi)
      {
        const auto& curF = fd.forces[fi];
        res(i) = -std::pow(curF.force()[2], 2)*muS +
            std::pow(curF.force()[0], 2) + std::pow(curF.force()[1], 2);
        ++i;
      }
    }
  }

private:
  PGData<Type>* pgdata_;
};

} // namespace pg
