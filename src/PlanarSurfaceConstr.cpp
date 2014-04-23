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

// associated header
#include "PlanarSurfaceConstr.h"

// include
// PG
#include "PGData.h"
#include "FillSparse.h"


namespace pg
{


/*
*                             PlanarPositionContactConstr
*/


PlanarPositionContactConstr::PlanarPositionContactConstr(PGData* pgdata, int bodyId,
    const sva::PTransformd& targetFrame,
    const sva::PTransformd& surfaceFrame)
  : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), 1, "PlanarPositionContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
  , targetFrame_(targetFrame)
  , surfaceFrame_(surfaceFrame)
  , jac_(pgdata->multibody(), bodyId, surfaceFrame.translation())
  , dotCache_(1, jac_.dof())
{}


PlanarPositionContactConstr::~PlanarPositionContactConstr() throw()
{ }


void PlanarPositionContactConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);

  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  res(0) = (pos.translation() - targetFrame_.translation()).dot(targetFrame_.rotation().row(2));
}


void PlanarPositionContactConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.reserve(jac_.dof());

  const Eigen::MatrixXd& mat = jac_.jacobian(pgdata_->multibody(), pgdata_->mbc());
  dotCache_.noalias() = targetFrame_.rotation().row(2)*mat.block(3, 0, 3, mat.cols());
  fullJacobianSparse(pgdata_->mb(), jac_, dotCache_, jac, {0, pgdata_->qParamsBegin()});
}


/*
*                             PlanarOrientationContactConstr
*/


PlanarOrientationContactConstr::PlanarOrientationContactConstr(PGData* pgdata, int bodyId,
    const sva::PTransformd& targetFrame,
    const sva::PTransformd& surfaceFrame,
    int axis)
  : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), 1, "PlanarPositionContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
  , targetFrame_(targetFrame)
  , surfaceFrame_(surfaceFrame)
  , axis_(axis)
  , jac_(pgdata->multibody(), bodyId, surfaceFrame.translation())
  , dotCache_(1, jac_.dof())
{}


PlanarOrientationContactConstr::~PlanarOrientationContactConstr() throw()
{ }


void PlanarOrientationContactConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);

  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  res(0) = (pos.rotation().row(axis_)).dot(targetFrame_.rotation().row(axis_));
}


void PlanarOrientationContactConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.reserve(jac_.dof());

  const Eigen::MatrixXd& mat = jac_.vectorJacobian(pgdata_->multibody(),
      pgdata_->mbc(), surfaceFrame_.rotation().row(axis_).transpose());
  dotCache_.noalias() = targetFrame_.rotation().row(axis_)*mat.block(3, 0, 3, mat.cols());
  fullJacobianSparse(pgdata_->mb(), jac_, dotCache_, jac, {0, pgdata_->qParamsBegin()});
}


/*
 *                             PlanarInclusionConstr
 */


PlanarInclusionConstr::PlanarInclusionConstr(PGData* pgdata, int bodyId,
    const sva::PTransformd& targetFrame,
    const std::vector<Eigen::Vector2d>& targetPoints,
    const sva::PTransformd& surfaceFrame,
    const std::vector<Eigen::Vector2d>& surfacePoints)
  : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), int(surfacePoints.size()*targetPoints.size()), "PlanarInclusionContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
  , targetFrame_(targetFrame)
  , targetPoints_(targetPoints)
  , targetVecNorm_(targetPoints.size())
  , surfacePoints_(surfacePoints.size())
  , jac_(pgdata->multibody(), bodyId)
  , transJac_(6, jac_.dof())
  , tJac_(1, jac_.dof())
  , bJac_(1, jac_.dof())
  , sumJac_(1, jac_.dof())
  , fullJac_(outputSize(), jac_.dof())
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
    surfacePoints_[i] = p*surfaceFrame;
  }
}


PlanarInclusionConstr::~PlanarInclusionConstr() throw()
{ }


void PlanarInclusionConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);

  int resIndex = 0;
  for(const sva::PTransformd& sp: surfacePoints_)
  {
    sva::PTransformd pos = sp*pgdata_->mbc().bodyPosW[bodyIndex_];
    Eigen::Vector3d posTargCoord = pos.translation() - targetFrame_.translation();
    double T = targetFrame_.rotation().row(0).dot(posTargCoord);
    double B = targetFrame_.rotation().row(1).dot(posTargCoord);
    Eigen::Vector2d vec(T, B);
    for(std::size_t i = 0; i < targetPoints_.size(); ++i)
    {
      const auto& n = targetVecNorm_[i];
      const auto& p = targetPoints_[i];
      res(resIndex) = n.x()*(vec.x() - p.x()) + n.y()*(vec.y() - p.y());
      ++resIndex;
    }
  }
}


void PlanarInclusionConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.reserve(outputSize()*jac_.dof());

  const Eigen::MatrixXd& jacMat = jac_.jacobian(pgdata_->multibody(), pgdata_->mbc());
  int resIndex = 0;
  for(const sva::PTransformd& sp: surfacePoints_)
  {
    jac_.translateJacobian(jacMat, pgdata_->mbc(), sp.translation(), transJac_);
    tJac_.noalias() = targetFrame_.rotation().row(0)*transJac_.block(3, 0, 3, jac_.dof());
    bJac_.noalias() = targetFrame_.rotation().row(1)*transJac_.block(3, 0, 3, jac_.dof());
    for(std::size_t i = 0; i < targetPoints_.size(); ++i)
    {
      const auto& n = targetVecNorm_[i];
      sumJac_.noalias() = n.x()*tJac_;
      sumJac_.noalias() += n.y()*bJac_;
      fullJac_.row(resIndex).noalias() = sumJac_;
      ++resIndex;
    }
  }
  fullJacobianSparse(pgdata_->mb(), jac_, fullJac_, jac, {0, pgdata_->qParamsBegin()});
}

} // pg
