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
#include "PositiveForceConstr.h"

// include
// PG
#include "PGData.h"

namespace pg
{


PositiveForceConstr::PositiveForceConstr(PGData* pgdata)
  : roboptim::DifferentiableFunction(pgdata->pbSize(), pgdata->nrForcePoints(), "PositiveForce")
  , pgdata_(pgdata)
  , jacPoints_(pgdata->nrForcePoints())
  , jacPointsMatTmp_(pgdata->nrForcePoints())
  , jacPointMatFull_(1, pgdata->mb().nrParams())
{
  std::size_t index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    for(std::size_t i = 0; i < fd.forces.size(); ++i)
    {
      jacPoints_[index] = rbd::Jacobian(pgdata_->mb(), fd.bodyId, fd.points[i].translation());
      jacPointsMatTmp_[index].resize(1, jacPoints_[index].dof());
      ++index;
    }
  }
}


PositiveForceConstr::~PositiveForceConstr() throw()
{ }


void PositiveForceConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);

  int index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    const sva::PTransformd& X_0_b = pgdata_->mbc().bodyPosW[fd.bodyIndex];
    for(std::size_t i = 0; i < fd.points.size(); ++i)
    {
      sva::PTransformd X_0_pi = fd.points[i]*X_0_b;
      res(index) = X_0_pi.rotation().row(2).dot(fd.forces[i].force());
      ++index;
    }
  }
}


void PositiveForceConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.setZero();

  int index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    const sva::PTransformd& X_0_b = pgdata_->mbc().bodyPosW[fd.bodyIndex];
    for(std::size_t i = 0; i < fd.points.size(); ++i)
    {
      sva::PTransformd X_0_pi = fd.points[i]*X_0_b;
      const Eigen::MatrixXd& jacPointMat =
          jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                           pgdata_->mbc(),
                                           fd.points[i].rotation().row(2).transpose())\
          .block(3, 0, 3, jacPoints_[index].dof());
      jacPointsMatTmp_[index].noalias() = fd.forces[i].force().transpose()*jacPointMat;
      jacPoints_[index].fullJacobian(pgdata_->mb(), jacPointsMatTmp_[index], jacPointMatFull_);
      jac.block(index, 0, 1, pgdata_->mb().nrParams()).noalias() = jacPointMatFull_;


      jac.block<1,3>(index, pgdata_->forceParamsBegin() + index*3).noalias() =
          X_0_pi.rotation().row(2);

      ++index;
    }
  }
}

} // pg

