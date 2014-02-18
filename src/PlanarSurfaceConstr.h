
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
    : parent_t(pgdata, pgdata->pbSize(), 1, "PlanarPositionContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , targetFrame_(targetFrame.cast<scalar_t>())
    , surfaceFrame_(surfaceFrame.cast<scalar_t>())
  {}
  ~PlanarPositionContactConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
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
class PlanarOrientationContactConstr : public AutoDiffFunction<Type, 3>
{
public:
  typedef AutoDiffFunction<Type, 3> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarOrientationContactConstr(PGData<Type>* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const sva::PTransformd& surfaceFrame,
      int axis)
    : parent_t(pgdata, pgdata->pbSize(), 3, "PlanarPositionContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , targetFrame_(targetFrame.cast<scalar_t>())
    , surfaceFrame_(surfaceFrame.cast<scalar_t>())
    , axis_(axis)
  {}
  ~PlanarOrientationContactConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
    sva::PTransform<scalar_t> pos = surfaceFrame_*pgdata_->fk().bodyPosW()[bodyIndex_];

    res(0) = (pos.rotation().row((axis_+1)%3)).dot(targetFrame_.rotation().row(axis_));
    res(1) = (pos.rotation().row((axis_+2)%3)).dot(targetFrame_.rotation().row(axis_));
    res(2) = (pos.rotation().row((axis_+3)%3)).dot(targetFrame_.rotation().row(axis_));
  }

private:
  PGData<Type>* pgdata_;

  int bodyIndex_;
  sva::PTransform<scalar_t> targetFrame_;
  sva::PTransform<scalar_t> surfaceFrame_;
  int axis_;
};



template<typename Type>
class PlanarInclusionConstr : public AutoDiffFunction<Type, Eigen::Dynamic>
{
public:
  typedef AutoDiffFunction<Type, Eigen::Dynamic> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarInclusionConstr(PGData<Type>* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const std::vector<Eigen::Vector2d>& targetPoints,
      const sva::PTransformd& surfaceFrame,
      const std::vector<Eigen::Vector2d>& surfacePoints)
    : parent_t(pgdata, pgdata->pbSize(), int(surfacePoints.size()*targetPoints.size()), "PlanarInclusionContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , targetFrame_(targetFrame.cast<scalar_t>())
    , targetPoints_(targetPoints)
    , targetVecNorm_(targetPoints.size())
    , surfacePoints_(surfacePoints.size())
  {
    assert(targetPoints.size() > 2);
    for(std::size_t i = 0; i < targetPoints.size(); ++i)
    {
      const Eigen::Vector2d& p1 = targetPoints[i];
      const Eigen::Vector2d& p2 = targetPoints[(i + 1) % targetPoints.size()];
      Eigen::Vector2d vec = p2 - p1;
      targetVecNorm_[i] = Eigen::Vector2d(-vec.y(), vec.x()).normalized();
    }

    for(std::size_t i = 0; i < surfacePoints.size(); ++i)
    {
      sva::PTransformd p(Eigen::Vector3d(surfacePoints[i].x(), surfacePoints[i].y(), 0.));
      surfacePoints_[i] = (p*surfaceFrame).cast<scalar_t>();
    }
  }
  ~PlanarInclusionConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
    int resIndex = 0;
    for(const sva::PTransform<scalar_t>& sp: surfacePoints_)
    {
      sva::PTransform<scalar_t> pos = sp*pgdata_->fk().bodyPosW()[bodyIndex_];
      Eigen::Vector3<scalar_t> posTargCoord = pos.translation() - targetFrame_.translation();
      scalar_t T = targetFrame_.rotation().row(0).dot(posTargCoord);
      scalar_t B = targetFrame_.rotation().row(1).dot(posTargCoord);
      Eigen::Matrix<scalar_t, 2, 1> vec(T, B);
      for(std::size_t i = 0; i < targetPoints_.size(); ++i)
      {
        const auto& n = targetVecNorm_[i];
        const auto& p = targetPoints_[i];
        res(resIndex) = n.x()*(vec.x() - p.x()) + n.y()*(vec.y() - p.y());
        ++resIndex;
      }
    }
  }

private:
  PGData<Type>* pgdata_;

  int bodyIndex_;
  sva::PTransform<scalar_t> targetFrame_;
  std::vector<Eigen::Vector2d> targetPoints_;
  std::vector<Eigen::Vector2d> targetVecNorm_;
  std::vector<sva::PTransform<scalar_t>> surfacePoints_; ///< Surface points in body coord.
};

} // namespace pg
