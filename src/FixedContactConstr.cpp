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
#include "FixedContactConstr.h"

// include
// PG
#include "PGData.h"

namespace pg
{


/*
 *                 FixedPositionContactConstr
 */


FixedPositionContactConstr::FixedPositionContactConstr(PGData* pgdata, int bodyId,
    const Eigen::Vector3d& target,
    const sva::PTransformd& surfaceFrame)
  : roboptim::DifferentiableFunction(pgdata->pbSize(), 3, "FixedPositionContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
  , target_(target)
  , surfaceFrame_(surfaceFrame)
  , jac_(pgdata->multibody(), bodyId, surfaceFrame.translation())
{}


FixedPositionContactConstr::~FixedPositionContactConstr() throw()
{ }


void FixedPositionContactConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);
  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  res = pos.translation() - target_;
}


void FixedPositionContactConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.setZero();
  const Eigen::MatrixXd& jacMat = jac_.jacobian(pgdata_->multibody(), pgdata_->mbc());
  jac_.fullJacobian(pgdata_->multibody(), jacMat.block(3, 0, 3, jacMat.cols()), jac);
}


/*
 *                      FixedOrientationContactConstr
 */


FixedOrientationContactConstr::FixedOrientationContactConstr(PGData* pgdata, int bodyId,
    const Eigen::Matrix3d& target,
    const sva::PTransformd& surfaceFrame)
  : roboptim::DifferentiableFunction(pgdata->pbSize(), 3, "FixedOrientationContact")
  , pgdata_(pgdata)
  , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
  , target_(target)
  , surfaceFrame_(surfaceFrame)
  , jac_(pgdata->multibody(), bodyId)
  , dotCache_(1, jac_.dof())
  , dotCacheFull_(1, pgdata_->multibody().nrDof())
{}


FixedOrientationContactConstr::~FixedOrientationContactConstr() throw()
{ }


void FixedOrientationContactConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);
  sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
  res(0) = pos.rotation().row(0).dot(target_.row(0));
  res(1) = pos.rotation().row(1).dot(target_.row(1));
  res(2) = pos.rotation().row(2).dot(target_.row(2));
}


template<typename Derived1, typename Derived2, typename Derived3>
void FixedOrientationContactConstr::dotDerivative(
    const Eigen::MatrixBase<Derived1>& posRow,
    const Eigen::MatrixBase<Derived2>& targetRow,
    Eigen::MatrixBase<Derived3> const & jac) const
{
  const Eigen::MatrixXd& mat =
    jac_.vectorJacobian(pgdata_->multibody(), pgdata_->mbc(), posRow.transpose());
  dotCache_.noalias() = targetRow*mat.block(3, 0, 3, mat.cols());
  jac_.fullJacobian(pgdata_->multibody(), dotCache_, dotCacheFull_);
  const_cast< Eigen::MatrixBase<Derived3>&>(jac)\
      .block(0, 0, 1, pgdata_->mb().nrParams()).noalias() = dotCacheFull_;
}


void FixedOrientationContactConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.setZero();

  dotDerivative(surfaceFrame_.rotation().row(0), target_.row(0), jac.row(0));
  dotDerivative(surfaceFrame_.rotation().row(1), target_.row(1), jac.row(1));
  dotDerivative(surfaceFrame_.rotation().row(2), target_.row(2), jac.row(2));
}


} // pg

