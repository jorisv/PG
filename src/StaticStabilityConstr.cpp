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
#include "StaticStabilityConstr.h"

// include
// PG
#include "PGData.h"


namespace pg
{


StaticStabilityConstr::StaticStabilityConstr(PGData* pgdata)
  : roboptim::DifferentiableFunction(pgdata->pbSize(), 6, "StaticStability")
  , pgdata_(pgdata)
  , gravityForce_(-pgdata->robotMass()*pgdata->gravity())
  , jacPoints_(pgdata->nrForcePoints())
  , jacFullMat_(3, pgdata->mb().nrParams())
  , T_com_fi_jac_(3, pgdata->mb().nrParams())
  , couple_jac_(3, pgdata->mb().nrParams())
{
  std::size_t index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    for(std::size_t i = 0; i < fd.forces.size(); ++i)
    {
      jacPoints_[index] = rbd::Jacobian(pgdata_->mb(), fd.bodyId, fd.points[i].translation());
      ++index;
    }
  }
}


StaticStabilityConstr::~StaticStabilityConstr() throw()
{ }


void StaticStabilityConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata_->x(x);

  res.head<3>().setZero();
  res.tail<3>() = gravityForce_;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    const sva::PTransformd& X_0_b = pgdata_->mbc().bodyPosW[fd.bodyIndex];
    for(std::size_t i = 0; i < fd.forces.size(); ++i)
    {
      sva::PTransformd X_0_pi = fd.points[i]*X_0_b;
      Eigen::Vector3d T_com_fi(X_0_pi.translation() - pgdata_->com());
      Eigen::Vector3d fi_world(fd.forces[i].force());
      res.head<3>() += T_com_fi.cross(fi_world);
      res.tail<3>() += fi_world;
    }
  }
}


void StaticStabilityConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  pgdata_->x(x);
  jac.setZero();
  std::size_t index = 0;
  for(const PGData::ForceData& fd: pgdata_->forceDatas())
  {
    const sva::PTransformd& X_0_b = pgdata_->mbc().bodyPosW[fd.bodyIndex];
    for(std::size_t i = 0; i < fd.forces.size(); ++i)
    {
      sva::PTransformd X_0_pi = fd.points[i]*X_0_b;
      Eigen::Vector3d T_com_fi(X_0_pi.translation() - pgdata_->com());

      // kinematic jacobian
      const Eigen::MatrixXd& jacP = jacPoints_[index].jacobian(pgdata_->mb(),
                                                                pgdata_->mbc());

      jacPoints_[index].fullJacobian(pgdata_->mb(), jacP.block(3, 0, 3, jacP.cols()), jacFullMat_);
      T_com_fi_jac_.noalias() = jacFullMat_ - pgdata_->comJac();
      couple_jac_.noalias() = sva::vector3ToCrossMatrix((-fd.forces[i].force()).eval())*T_com_fi_jac_;

      // couple
      jac.block(0, 0, 3, pgdata_->mb().nrParams()).noalias() += couple_jac_;
      // force
      // Zero


      // force jacobian
      // couple
      jac.block<3,3>(0, pgdata_->forceParamsBegin() + index*3).noalias() =
          sva::vector3ToCrossMatrix(T_com_fi);
      // force
      jac.block<3,3>(3, pgdata_->forceParamsBegin() + index*3).noalias() =
          Eigen::Matrix3d::Identity();
      ++index;
    }
  }
}

} // pg
