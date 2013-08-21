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
class FixedContactConstr : public AutoDiffFunction<Type, 3>
{
public:
  typedef AutoDiffFunction<Type, 3> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  FixedContactConstr(PGData<Type>* pgdata, int bodyId, const Eigen::Vector3d& pos)
    : parent_t(pgdata->pbSize(), 3, "FixedContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , pos_(pos.cast<scalar_t>())
  {}
  ~FixedContactConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& x) const
  {
    pgdata_->x(x);
    res = pgdata_->fk().bodyPosW()[bodyIndex_].translation() - pos_;
  }

private:
  PGData<Type>* pgdata_;

  int bodyIndex_;
  Eigen::Vector3<scalar_t> pos_;
};

} // namespace pg
