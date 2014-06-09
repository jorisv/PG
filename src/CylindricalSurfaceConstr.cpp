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
#include "CylindricalSurfaceConstr.h"

// include
// PG
#include "PGData.h"
#include "FillSparse.h"


namespace pg
{


/*
 *                        CylindricalPositionConstr
 */


CylindricalPositionConstr::CylindricalPositionConstr(
  PGData* pgdata, int bodyId,
  const sva::PTransformd& targetFrame,
  const sva::PTransformd& surfaceFrame)
  : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), 3, "FreeGripperPositionContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
  , targetFrame_(targetFrame)
  , surfaceFrame_(surfaceFrame)
  , jac_(pgdata->mb(), bodyId, surfaceFrame.translation())
  , jacMat_(3, jac_.dof())
{}


CylindricalPositionConstr::~CylindricalPositionConstr()
{}


void CylindricalPositionConstr::impl_compute(result_t& res, const argument_t& x) const
{
  pgdata_->x(x);

  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  Eigen::Vector3d vec = pos.translation() - targetFrame_.translation();

  res(0) = targetFrame_.rotation().row(0).dot(vec);
  res(1) = targetFrame_.rotation().row(1).dot(vec);
  res(2) = targetFrame_.rotation().row(2).dot(vec);
}


void CylindricalPositionConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const
{
  jac.reserve(3*jac_.dof());

  pgdata_->x(x);

  const Eigen::MatrixXd& jacMat = jac_.jacobian(pgdata_->mb(), pgdata_->mbc());

  jacMat_.row(0).noalias() = targetFrame_.rotation().row(0)*jacMat.block(3, 0, 3, jac_.dof());
  jacMat_.row(1).noalias() = targetFrame_.rotation().row(1)*jacMat.block(3, 0, 3, jac_.dof());
  jacMat_.row(2).noalias() = targetFrame_.rotation().row(2)*jacMat.block(3, 0, 3, jac_.dof());
  fullJacobianSparse(pgdata_->mb(), jac_, jacMat_,
                     jac, {0, pgdata_->qParamsBegin()});
}


/*
 *                        CylindricalNVecConstr
 */


CylindricalNVecConstr::CylindricalNVecConstr(
  PGData* pgdata, int bodyId,
  const sva::PTransformd& targetFrame,
  const sva::PTransformd& surfaceFrame)
  : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), 1, "FreeGripperNVecContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->mb().bodyIndexById(bodyId))
  , targetFrame_(targetFrame)
  , surfaceFrame_(surfaceFrame)
  , jac_(pgdata->mb(), bodyId)
  , jacMat_(1, jac_.dof())
{}


CylindricalNVecConstr::~CylindricalNVecConstr()
{}


void CylindricalNVecConstr::impl_compute(result_t& res, const argument_t& x) const
{
  pgdata_->x(x);

  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  Eigen::Vector3d vec = (targetFrame_.translation() - pos.translation());
  double dot = vec.dot(pos.rotation().row(2));
  double sign = std::copysign(1., dot);
  res(0) = sign*std::pow(dot, 2) - vec.squaredNorm();
}


void CylindricalNVecConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const
{
  jac.reserve(jac_.dof());

  pgdata_->x(x);

  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  Eigen::Vector3d vec = (targetFrame_.translation() - pos.translation());
  double dot = vec.dot(pos.rotation().row(2));
  double sign = std::copysign(1., dot);
  double vecNDot = sign*2.*vec.dot(pos.rotation().row(2));

  const Eigen::MatrixXd& jacVecNMat =
      jac_.vectorJacobian(pgdata_->mb(), pgdata_->mbc(),
                          surfaceFrame_.rotation().row(2).transpose());

  jacMat_.row(0).noalias() = vecNDot*vec.transpose()*jacVecNMat.block(3, 0, 3, jac_.dof());

  jac_.point(surfaceFrame_.translation());
  const Eigen::MatrixXd& jacMat = jac_.jacobian(pgdata_->mb(), pgdata_->mbc());
  jac_.point(Eigen::Vector3d::Zero());

  jacMat_.row(0).noalias() -= vecNDot*pos.rotation().row(2)*jacMat.block(3, 0, 3, jac_.dof());
  jacMat_.row(0).noalias() += (2.*vec.transpose())*jacMat.block(3, 0, 3, jac_.dof());

  fullJacobianSparse(pgdata_->mb(), jac_, jacMat_,
                     jac, {0, pgdata_->qParamsBegin()});
}


} // pg
