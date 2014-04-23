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
#include "StdCostFunc.h"

// include
// boost
#include <boost/math/constants/constants.hpp>

// PG
#include "PGData.h"
#include "ConfigStruct.h"
#include "FillSparse.h"


namespace pg
{


StdCostFunc::StdCostFunc(std::vector<PGData>& pgdatas,
                         const std::vector<RobotConfig>& robotConfigs,
                         const std::vector<RunConfig>& runConfigs)
  : roboptim::DifferentiableSparseFunction(pgdatas[0].pbSize(), 1, "StdCostFunc")
  , robotDatas_(pgdatas.size())
{
  for(std::size_t robotIndex = 0; robotIndex < pgdatas.size(); ++robotIndex)
  {
    const RobotConfig& robotConfig = robotConfigs[robotIndex];
    const RunConfig& runConfig = runConfigs[robotIndex];
    PGData& pgdata = pgdatas[robotIndex];

    RobotData data;
    data.pgdata = &pgdata;
    data.tq = runConfig.targetQ;
    data.postureScale = robotConfig.postureScale;
    data.torqueScale = robotConfig.torqueScale;
    data.forceScale = robotConfig.forceScale;
    data.ellipseScale = robotConfig.ellipseCostScale;

    for(std::size_t i = 0; i < robotConfig.bodyPosTargets.size(); ++i)
    {
      const BodyPositionTarget& bodyPosTarget = robotConfig.bodyPosTargets[i];
      rbd::Jacobian jac(pgdata.mb(), bodyPosTarget.bodyId);
      Eigen::MatrixXd jacMat(1, jac.dof());
      Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull(1, pgdata.pbSize());
      jacMatFull.reserve(jac.dof());
      data.bodyPosTargets.push_back({pgdata.mb().bodyIndexById(bodyPosTarget.bodyId),
                                     bodyPosTarget.target,
                                     bodyPosTarget.scale,
                                     jac, jacMat, jacMatFull});
    }

    for(std::size_t i = 0; i < robotConfig.bodyOriTargets.size(); ++i)
    {
      const BodyOrientationTarget& bodyOriTarget = robotConfig.bodyOriTargets[i];
      rbd::Jacobian jac(pgdata.mb(), bodyOriTarget.bodyId);
      Eigen::MatrixXd jacMat(1, jac.dof());
      Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull(1, pgdata.pbSize());
      jacMatFull.reserve(jac.dof());
      data.bodyOriTargets.push_back({pgdata.mb().bodyIndexById(bodyOriTarget.bodyId),
                                     bodyOriTarget.target,
                                     bodyOriTarget.scale,
                                     jac, jacMat, jacMatFull});
    }

    for(std::size_t i = 0; i < robotConfig.forceContactsMin.size(); ++i)
    {
      std::size_t gradientPos = pgdata.forceParamsBegin();
      for(std::size_t j = 0; j < robotConfig.forceContacts.size(); ++j)
      {
        // we don't break since it could be many contact on the same body
        if(robotConfig.forceContactsMin[i].bodyId == robotConfig.forceContacts[j].bodyId)
        {
          data.forceContactsMin.push_back({j, gradientPos,
                                           robotConfig.forceContactsMin[i].scale});
        }
        gradientPos += robotConfig.forceContacts[j].points.size()*3;
      }
    }

    robotDatas_[robotIndex] = std::move(data);
  }
}


void StdCostFunc::impl_compute(result_t& res, const argument_t& x) const throw()
{
  res(0) = 0.;

  for(RobotData& rd: robotDatas_)
  {
    rd.pgdata->x(x);

    // compute posture task
    double posture = 0.;
    if(rd.postureScale > 0.)
    {
      const std::vector<std::vector<double>>& q = rd.pgdata->mbc().q;
      for(int i = 0; i < rd.pgdata->mb().nrJoints(); ++i)
      {
        if(rd.pgdata->mb().joint(i).params() == 1)
        {
          posture += std::pow(q[i][0] - rd.tq[i][0], 2);
        }
      }
    }

    // compute torque task
    /*
    scalar_t torque = scalar_t(0., Eigen::VectorXd::Zero(this->inputSize()));
    if(torqueScale_ > 0.)
    {
      const auto& torqueVec = rd.pgdata->id().torque();
      for(int i = 0; i < rd.pgdata->multibody().nrJoints(); ++i)
      {
        for(int j = 0; j < rd.pgdata->multibody().joint(i).dof(); ++j)
        {
          torque += std::pow(torqueVec[i](j), 2);
        }
      }
    }
    */

    // compute force task
    double force = 0.;
    if(rd.forceScale > 0.)
    {
      for(const auto& fd: rd.pgdata->forceDatas())
      {
        Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
        for(const sva::ForceVecd& fv: fd.forces)
        {
          forceTmp += fv.force();
        }
        force += forceTmp.squaredNorm()*rd.forceScale;
      }
    }

    double pos = 0.;
    for(const BodyPositionTargetData& bp: rd.bodyPosTargets)
    {
      pos += (rd.pgdata->mbc().bodyPosW[bp.bodyIndex].translation() - bp.target).squaredNorm()*
          bp.scale;
    }

    double ori = 0.;
    for(const BodyOrientationTargetData& bo: rd.bodyOriTargets)
    {
      ori += sva::rotationError(bo.target,
        rd.pgdata->mbc().bodyPosW[bo.bodyIndex].rotation(), 1e-7).squaredNorm()*
          bo.scale;
    }

    double forceMin = 0.;
    for(const ForceContactMinimizationData& fcmd: rd.forceContactsMin)
    {
      const auto& forceData = rd.pgdata->forceDatas()[fcmd.forcePos];
      Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
      for(const sva::ForceVecd& fv: forceData.forces)
      {
        forceTmp += fv.force();
      }
      forceMin += forceTmp.squaredNorm()*fcmd.scale;
    }

    //Compute ellipse contact cost
    /*
    scalar_t ellipses = scalar_t(0., Eigen::VectorXd::Zero(this->inputSize()));
    if(ellipseScale_ > 0.)
    {
      for(std::size_t i = 0; i < rd.pgdata->ellipseDatas().size(); ++i)
      {
        ellipses += -boost::math::constants::pi<double>()*rd.pgdata->ellipseDatas()[i].r1 * rd.pgdata->ellipseDatas()[i].r2;
        std::cout << "ellipses cost: " << ellipses << std::endl;
      }
    }
    */

    res(0) += posture*rd.postureScale + force +
        pos + ori + forceMin;
  }
}

void StdCostFunc::impl_gradient(gradient_t& gradient,
    const argument_t& x, size_type /* functionId */) const throw()
{
  gradient.reserve(robotDatas_[0].pgdata->pbSize());
  gradient.setZero();

  for(RobotData& rd: robotDatas_)
  {
    rd.pgdata->x(x);

    if(rd.postureScale > 0.)
    {
      int index = rd.pgdata->qParamsBegin();
      const std::vector<std::vector<double>>& q = rd.pgdata->mbc().q;
      double coef = 2.*rd.postureScale;
      for(int i = 0; i < rd.pgdata->multibody().nrJoints(); ++i)
      {
        if(rd.pgdata->multibody().joint(i).params() == 1)
        {
          gradient.coeffRef(index) += coef*(q[i][0] - rd.tq[i][0]);
        }
        index += rd.pgdata->mb().joint(i).dof();
      }
    }


    if(rd.forceScale > 0.)
    {
      int index = rd.pgdata->forceParamsBegin();
      for(const auto& fd: rd.pgdata->forceDatas())
      {
        Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
        for(const sva::ForceVecd& fv: fd.forces)
        {
          forceTmp += fv.force();
        }
        forceTmp *= 2.*rd.forceScale;

        for(std::size_t i = 0; i < fd.forces.size(); ++i)
        {
          gradient.coeffRef(index + 0) += forceTmp.x();
          gradient.coeffRef(index + 1) += forceTmp.y();
          gradient.coeffRef(index + 2) += forceTmp.z();
          index += 3;
        }
      }
    }


    for(const ForceContactMinimizationData& fcmd: rd.forceContactsMin)
    {
      const auto& forceData = rd.pgdata->forceDatas()[fcmd.forcePos];
      Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
      int index = int(fcmd.gradientPos);

      for(const sva::ForceVecd& fv: forceData.forces)
      {
        forceTmp += fv.force();
      }
      forceTmp *= 2.*fcmd.scale;

      for(std::size_t i = 0; i < forceData.forces.size(); ++i)
      {
        gradient.coeffRef(index + 0) += forceTmp.x();
        gradient.coeffRef(index + 1) += forceTmp.y();
        gradient.coeffRef(index + 2) += forceTmp.z();
        index += 3;
      }
    }


    for(BodyPositionTargetData& bp: rd.bodyPosTargets)
    {
      Eigen::Vector3d error(rd.pgdata->mbc().bodyPosW[bp.bodyIndex].translation()
        - bp.target);

      const Eigen::MatrixXd& jacMat = bp.jac.jacobian(rd.pgdata->mb(), rd.pgdata->mbc());
      bp.jacMat.noalias() = (bp.scale*2.*error.transpose())*jacMat.block(3, 0, 3, bp.jac.dof());
      updateFullJacobianSparse(rd.pgdata->mb(), bp.jac, bp.jacMat, bp.jacMatFull,
                               {0, rd.pgdata->qParamsBegin()});
      gradient += bp.jacMatFull.transpose();
    }


    for(BodyOrientationTargetData& bo: rd.bodyOriTargets)
    {
      Eigen::Vector3d error(sva::rotationError(bo.target,
        rd.pgdata->mbc().bodyPosW[bo.bodyIndex].rotation(), 1e-7));

      const Eigen::MatrixXd& jacMat = bo.jac.jacobian(rd.pgdata->mb(), rd.pgdata->mbc());
      bo.jacMat.noalias() = (bo.scale*2.*error.transpose())*jacMat.block(0, 0, 3, bo.jac.dof());
      updateFullJacobianSparse(rd.pgdata->mb(), bo.jac, bo.jacMat, bo.jacMatFull,
                               {0, rd.pgdata->qParamsBegin()});
      gradient += bo.jacMatFull.transpose();
    }
  }
}


} // pg
