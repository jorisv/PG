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
#include "FrictionConeConstr.h"

// PG
#include "ConfigStruct.h"
#include "PGData.h"
#include "FillSparse.h"

namespace pg
{


FrictionConeConstr::FrictionConeConstr(PGData* pgdata)
  : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), pgdata->nrForcePoints(), "FrictionConeConstr")
  , pgdata_(pgdata)
  , nrNonZero_(0)
  , jacPoints_(pgdata->nrForcePoints())
  , jacPointsMatTmp_(pgdata->nrForcePoints())
{
  std::size_t index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    for(std::size_t i = 0; i < fd.forces.size(); ++i)
    {
      jacPoints_[index] = rbd::Jacobian(pgdata_->mb(), fd.bodyId, fd.points[i].translation());
      jacPointsMatTmp_[index].resize(1, jacPoints_[index].dof());
      // jacobian dof + force vec
      nrNonZero_ += jacPoints_[index].dof() + 3;
      ++index;
    }
  }
}


FrictionConeConstr::~FrictionConeConstr() throw()
{ }


void FrictionConeConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);

  int index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    const sva::PTransformd& X_0_b = pgdata_->mbc().bodyPosW[fd.bodyIndex];
    for(std::size_t i = 0; i < fd.points.size(); ++i)
    {
      sva::PTransformd X_0_pi = fd.points[i]*X_0_b;
      Eigen::Vector3d fBody(X_0_pi.rotation()*fd.forces[i].force());
      res(index) = -std::pow(fBody.z()*fd.mu, 2) +
          std::pow(fBody.x(), 2) + std::pow(fBody.y(), 2);
      ++index;
    }
  }
}


void FrictionConeConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.reserve(nrNonZero_);

  int index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    const sva::PTransformd& X_0_b = pgdata_->mbc().bodyPosW[fd.bodyIndex];
    double muSquare = std::pow(fd.mu, 2);
    for(std::size_t i = 0; i < fd.points.size(); ++i)
    {
      sva::PTransformd X_0_pi = fd.points[i]*X_0_b;
      Eigen::Vector3d fBody(X_0_pi.rotation()*fd.forces[i].force());

      // Z axis
      //                dq
      // -(mu*forceB.z())^2 => -2*mu^2*forceB.z()*forceW*jacZ
      //
      // X axis
      //               dq
      // forceB.x()**2 => 2*forceB.x()*forceW*jacX
      //
      // Y axis
      //               dq
      // forceB.y()**2 => 2*forceB.y()*forceW*jacY
      const Eigen::MatrixXd& jacPointVecZMat =
          jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                           pgdata_->mbc(),
                                           fd.points[i].rotation().row(2).transpose())\
          .block(3, 0, 3, jacPoints_[index].dof());
      jacPointsMatTmp_[index].row(0).noalias() =
          (-2.*muSquare*fBody.z())*(fd.forces[i].force().transpose()*jacPointVecZMat);

      const Eigen::MatrixXd& jacPointVecXMat =
          jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                           pgdata_->mbc(),
                                           fd.points[i].rotation().row(0).transpose())\
          .block(3, 0, 3, jacPoints_[index].dof());
      jacPointsMatTmp_[index].noalias() +=
          (2.*fBody.x())*(fd.forces[i].force().transpose()*jacPointVecXMat);

      const Eigen::MatrixXd& jacPointVecYMat =
          jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                           pgdata_->mbc(),
                                           fd.points[i].rotation().row(1).transpose())\
          .block(3, 0, 3, jacPoints_[index].dof());
      jacPointsMatTmp_[index].noalias() +=
          (2.*fBody.y())*(fd.forces[i].force().transpose()*jacPointVecYMat);

      fullJacobianSparse(pgdata_->mb(), jacPoints_[index], jacPointsMatTmp_[index],
        jac, {index, pgdata_->qParamsBegin()});

      // Z axis
      //                dforceW
      // -(mu*forceB.z())^2 => -2*mu^2*forceB.z()*X_0_pi_rowZ
      //
      // X axis
      //               dforceW
      // forceB.x()**2 => 2*forceB.x()*X_0_pi_rowX
      //
      // Y axis
      //               dforceW
      // forceB.y()**2 => 2*forceB.y()*X_0_pi_rowY
      Eigen::RowVector3d fDiff = (-2.*muSquare*fBody.z())*X_0_pi.rotation().row(2);
      fDiff.noalias() += (2.*fBody.x())*X_0_pi.rotation().row(0);
      fDiff.noalias() += (2.*fBody.y())*X_0_pi.rotation().row(1);
      int indexCols = pgdata_->forceParamsBegin() + index*3;
      jac.insert(index, indexCols + 0) = fDiff(0);
      jac.insert(index, indexCols + 1) = fDiff(1);
      jac.insert(index, indexCols + 2) = fDiff(2);

      ++index;
    }
  }
}


} // pg
