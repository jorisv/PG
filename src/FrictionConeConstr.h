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
#include "PGData.h"

namespace pg
{

class FrictionConeConstr : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  FrictionConeConstr(PGData* pgdata)
    : roboptim::DifferentiableFunction(pgdata->pbSize(), pgdata->nrForcePoints(), "FrictionConeConstr")
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
  ~FrictionConeConstr() throw()
  { }


  void impl_compute(result_t& res, const argument_t& x) const throw()
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
        res(index) = -fBody.z()*fd.mu +
            std::sqrt(std::pow(fBody.x(), 2) + std::pow(fBody.y(), 2));
        ++index;
      }
    }
  }


  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
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
        Eigen::Vector3d fBody(X_0_pi.rotation()*fd.forces[i].force());

        // Z axis
        //                dq
        // -mu*forceB.z() => -mu*forceW*jacZ
        const Eigen::MatrixXd& jacPointVecZMat =
            jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                             pgdata_->mbc(),
                                             fd.points[i].rotation().row(2).transpose())\
            .block(3, 0, 3, jacPoints_[index].dof());
        jacPointsMatTmp_[index].row(0).noalias() = -fd.mu*fd.forces[i].force().transpose()*jacPointVecZMat;
        jacPoints_[index].fullJacobian(pgdata_->mb(), jacPointsMatTmp_[index], jacPointMatFull_);
        jac.block(index, 0, 1, pgdata_->mb().nrParams()).noalias() = jacPointMatFull_;

        // X axis
        //               dq
        // forceB.x()**2 => 2*forceB.x()*forceW*jacX
        //
        // Y axis
        //               dq
        // forceB.y()**2 => 2*forceB.y()*forceW*jacY
        //
        // sqrt(forceB.x()**2 + forceB.y()**2)
        // dq =>
        // d forceB.x()**2   d forceB.y()**2
        // --------------- + ---------------
        //     dq                dq
        // ---------------------------------
        // 2*sqrt(forceB.x()**2 + forceB.y()**2)

        double xyNorm = std::sqrt(std::pow(fBody.x(), 2) + std::pow(fBody.y(), 2));
        const Eigen::MatrixXd& jacPointVecXMat =
            jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                             pgdata_->mbc(),
                                             fd.points[i].rotation().row(0).transpose())\
            .block(3, 0, 3, jacPoints_[index].dof());
        jacPointsMatTmp_[index].noalias() =
            (fBody.x()/xyNorm)*fd.forces[i].force().transpose()*jacPointVecXMat;

        const Eigen::MatrixXd& jacPointVecYMat =
            jacPoints_[index].vectorJacobian(pgdata_->mb(),
                                             pgdata_->mbc(),
                                             fd.points[i].rotation().row(1).transpose())\
            .block(3, 0, 3, jacPoints_[index].dof());
        jacPointsMatTmp_[index].noalias() +=
            (fBody.y()/xyNorm)*fd.forces[i].force().transpose()*jacPointVecYMat;

        jacPoints_[index].fullJacobian(pgdata_->mb(), jacPointsMatTmp_[index], jacPointMatFull_);
        jac.block(index, 0, 1, pgdata_->mb().nrParams()).noalias() += jacPointMatFull_;

        // Z axis
        //                dforceW
        // -mu*forceB.z() => -mu*X_0_pi_rowZ
        jac.block<1,3>(index, pgdata_->forceParamsBegin() + index*3).noalias() =
            -X_0_pi.rotation().row(2)*fd.mu;

        // X axis
        //               dforceW
        // forceB.x()**2 => 2*forceB.x()*X_0_pi_rowX
        //
        // Y axis
        //               dforceW
        // forceB.y()**2 => 2*forceB.y()*X_0_pi_rowY
        //
        // sqrt(forceB.x()**2 + forceB.y()**2)
        // dforceW =>
        // d forceB.x()**2   d forceB.y()**2
        // --------------- + ---------------
        //     dforceW                dforceW
        // ---------------------------------
        // 2*sqrt(forceB.x()**2 + forceB.y()**2)
        jac.block<1,3>(index, pgdata_->forceParamsBegin() + index*3).noalias() +=
            (fBody.x()/xyNorm)*X_0_pi.rotation().row(0);
        jac.block<1,3>(index, pgdata_->forceParamsBegin() + index*3).noalias() +=
            (fBody.y()/xyNorm)*X_0_pi.rotation().row(1);

        ++index;
      }
    }
  }

  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;

  mutable std::vector<rbd::Jacobian> jacPoints_;
  mutable std::vector<Eigen::MatrixXd> jacPointsMatTmp_;
  mutable Eigen::MatrixXd jacPointMatFull_;
};

} // namespace pg
