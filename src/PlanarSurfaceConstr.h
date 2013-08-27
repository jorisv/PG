
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
class PlanarPositionContactConstr : public AutoDiffFunction<Type, 1>
{
public:
  typedef AutoDiffFunction<Type, 1> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarPositionContactConstr(PGData<Type>* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const sva::PTransformd& surfaceFrame)
    : parent_t(pgdata->pbSize(), 1, "PlanarPositionContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , targetFrame_(targetFrame.cast<scalar_t>())
    , surfaceFrame_(surfaceFrame.cast<scalar_t>())
  {}
  ~PlanarPositionContactConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& x) const
  {
    pgdata_->x(x);
    sva::PTransform<scalar_t> pos = surfaceFrame_*pgdata_->fk().bodyPosW()[bodyIndex_];

    res(0) = (pos.translation() - targetFrame_.translation()).dot(targetFrame_.rotation().row(2));
  }

private:
  PGData<Type>* pgdata_;

  int bodyIndex_;
  sva::PTransform<scalar_t> targetFrame_;
  sva::PTransform<scalar_t> surfaceFrame_;
};



template<typename Type>
class PlanarOrientationContactConstr : public AutoDiffFunction<Type, 1>
{
public:
  typedef AutoDiffFunction<Type, 1> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarOrientationContactConstr(PGData<Type>* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const sva::PTransformd& surfaceFrame)
    : parent_t(pgdata->pbSize(), 1, "PlanarPositionContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , targetFrame_(targetFrame.cast<scalar_t>())
    , surfaceFrame_(surfaceFrame.cast<scalar_t>())
  {}
  ~PlanarOrientationContactConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& x) const
  {
    pgdata_->x(x);
    sva::PTransform<scalar_t> pos = surfaceFrame_*pgdata_->fk().bodyPosW()[bodyIndex_];

    res(0) = (pos.rotation().row(2)).dot(targetFrame_.rotation().row(2));
  }

private:
  PGData<Type>* pgdata_;

  int bodyIndex_;
  sva::PTransform<scalar_t> targetFrame_;
  sva::PTransform<scalar_t> surfaceFrame_;
};

} // namespace pg
