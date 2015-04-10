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
  , scale_(1./double(robotConfigs.size()))
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
      for(std::size_t j = 0; j < pgdata.forceDatas().size(); ++j)
      {
        const PGData::ForceData& forceData = pgdata.forceDatas()[j];
        // we don't break since it could be many contact on the same body
        if(robotConfig.forceContactsMin[i].bodyId == forceData.bodyId)
        {
          data.forceContactsMin.push_back({j, gradientPos,
                                           robotConfig.forceContactsMin[i].scale});
        }
        gradientPos += forceData.points.size()*3;
      }
    }

    for(std::size_t i = 0; i < robotConfig.torqueContactsMin.size(); ++i)
    {
      std::size_t gradientPos = pgdata.forceParamsBegin();
      const TorqueContactMinimization& tcm = robotConfig.torqueContactsMin[i];
      int bodyIndex = pgdata.mb().bodyIndexById(tcm.bodyId);
      for(std::size_t j = 0; j < pgdata.forceDatas().size(); ++j)
      {
        const PGData::ForceData& forceData = pgdata.forceDatas()[j];
        // we don't break since it could be many contact on the same body
        if(tcm.bodyId == forceData.bodyId)
        {
          std::vector<Eigen::Vector3d> levers;
          levers.reserve(forceData.points.size());
          for(const sva::PTransformd& p: forceData.points)
          {
            levers.push_back(p.translation() - tcm.origin);
          }
          rbd::Jacobian jac(pgdata.mb(), tcm.bodyId);
          Eigen::MatrixXd jacMat(1, jac.dof());
          Eigen::MatrixXd jacMatTmp(3, jac.dof());
          Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull(1, pgdata.pbSize());
          jacMatFull.reserve(jac.dof());
          data.torqueContactsMin.push_back({bodyIndex, forceData.points, std::move(levers),
                                            tcm.axis, jac, jacMat, jacMatTmp, jacMatFull,
                                            j, gradientPos, tcm.scale});
        }
        gradientPos += forceData.points.size()*3;
      }
    }

    for(std::size_t i = 0; i < robotConfig.normalForceTargets.size(); ++i)
    {
      std::size_t gradientPos = pgdata.forceParamsBegin();
      for(std::size_t j = 0; j < pgdata.forceDatas().size(); ++j)
      {
        const PGData::ForceData& forceData = pgdata.forceDatas()[j];
        // we don't break since it could be many contact on the same body
        if(robotConfig.normalForceTargets[i].bodyId == forceData.bodyId)
        {
          rbd::Jacobian jac(pgdata.mb(), robotConfig.normalForceTargets[i].bodyId);
          Eigen::MatrixXd jacMat(1, jac.dof());
          Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull(1, pgdata.pbSize());
          jacMatFull.reserve(jac.dof());
          data.normalForceTargets.push_back({j, gradientPos,
                                             robotConfig.normalForceTargets[i].target,
                                             jac, jacMat, jacMatFull,
                                             robotConfig.normalForceTargets[i].scale});
        }
        gradientPos += forceData.points.size()*3;
      }
    }

    for(std::size_t i = 0; i < robotConfig.tanForceMin.size(); ++i)
    {
      std::size_t gradientPos = pgdata.forceParamsBegin();
      for(std::size_t j = 0; j < pgdata.forceDatas().size(); ++j)
      {
        const PGData::ForceData& forceData = pgdata.forceDatas()[j];
        // we don't break since it could be many contact on the same body
        if(robotConfig.tanForceMin[i].bodyId == forceData.bodyId)
        {
          rbd::Jacobian jac(pgdata.mb(), robotConfig.tanForceMin[i].bodyId);
          Eigen::MatrixXd jacMat(1, jac.dof());
          Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull(1, pgdata.pbSize());
          jacMatFull.reserve(jac.dof());
          data.tanForceMin.push_back({j, gradientPos,
                                      jac, jacMat, jacMatFull,
                                      robotConfig.tanForceMin[i].scale});
        }
        gradientPos += forceData.points.size()*3;
      }
    }

    robotDatas_[robotIndex] = std::move(data);
  }
}


void StdCostFunc::impl_compute(result_t& res, const argument_t& x) const
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
        for(const sva::ForceVecd& fv: fd.forces)
        {
          force += fv.force().squaredNorm()*rd.forceScale;
        }
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
      for(const sva::ForceVecd& fv: forceData.forces)
      {
        forceMin += fv.force().squaredNorm()*fcmd.scale;
      }
    }

    double torqueContactMin = 0.;
    for(const TorqueContactMinimizationData& tcmd: rd.torqueContactsMin)
    {
      const auto& forceData = rd.pgdata->forceDatas()[tcmd.forcePos];
      const sva::PTransformd& X_0_b = rd.pgdata->mbc().bodyPosW[tcmd.bodyIndex];
      for(std::size_t pi = 0; pi < tcmd.points.size(); ++pi)
      {
        Eigen::Vector3d fBody(X_0_b.rotation()*forceData.forces[pi].force());
        torqueContactMin += std::pow(tcmd.axis.dot(tcmd.levers[pi].cross(fBody)), 2)*tcmd.scale;
      }
    }

    double normalForceTarget = 0.;
    for(const NormalForceTargetData& nft: rd.normalForceTargets)
    {
      const auto& forceData = rd.pgdata->forceDatas()[nft.forcePos];
      const sva::PTransformd& X_0_b = rd.pgdata->mbc().bodyPosW[forceData.bodyIndex];
      double nForceSum = 0.;
      for(std::size_t pi = 0; pi < forceData.points.size(); ++pi)
      {
        // force in the point pi frame
        sva::PTransformd X_0_pi = forceData.points[pi]*X_0_b;
        nForceSum += (X_0_pi.rotation()*forceData.forces[pi].force()).z();
      }
      normalForceTarget += std::pow(nForceSum - nft.target, 2)*nft.scale;
    }

    double tanForceMin = 0.;
    for(const TangentialForceMinimizationData& tfm: rd.tanForceMin)
    {
      const auto& forceData = rd.pgdata->forceDatas()[tfm.forcePos];
      const sva::PTransformd& X_0_b = rd.pgdata->mbc().bodyPosW[forceData.bodyIndex];
      for(std::size_t pi = 0; pi < forceData.points.size(); ++pi)
      {
        // force in the point pi frame
        sva::PTransformd X_0_pi = forceData.points[pi]*X_0_b;
        Eigen::Vector3d fPoint(X_0_pi.rotation()*forceData.forces[pi].force());
        tanForceMin += Eigen::Vector2d(fPoint.x(), fPoint.y()).squaredNorm()*tfm.scale;
      }
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

    res(0) += (posture*rd.postureScale + force +
        pos + ori + forceMin + torqueContactMin + normalForceTarget +
        tanForceMin)*scale_;
  }
}

void StdCostFunc::impl_gradient(gradient_t& gradient,
    const argument_t& x, size_type /* functionId */) const
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
      double coef = 2.*rd.postureScale*scale_;
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
        for(std::size_t i = 0; i < fd.forces.size(); ++i)
        {
          Eigen::Vector3d forceTmp = fd.forces[i].force()*rd.forceScale*scale_*2.;

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
      int index = int(fcmd.gradientPos);

      for(std::size_t i = 0; i < forceData.forces.size(); ++i)
      {
        Eigen::Vector3d forceTmp = forceData.forces[i].force()*fcmd.scale*scale_*2.;

        gradient.coeffRef(index + 0) += forceTmp.x();
        gradient.coeffRef(index + 1) += forceTmp.y();
        gradient.coeffRef(index + 2) += forceTmp.z();
        index += 3;
      }
    }


    for(TorqueContactMinimizationData& tcmd: rd.torqueContactsMin)
    {
      const auto& forceData = rd.pgdata->forceDatas()[tcmd.forcePos];
      const sva::PTransformd& X_0_b = rd.pgdata->mbc().bodyPosW[tcmd.bodyIndex];

      int forceIndex = int(tcmd.gradientPos);
      for(std::size_t pi = 0; pi < tcmd.points.size(); ++pi)
      {
        Eigen::Vector3d fBody(X_0_b.rotation()*forceData.forces[pi].force());
        Eigen::Vector3d lever = tcmd.levers[pi].cross(fBody);
        Eigen::Matrix3d crossMat = sva::vector3ToCrossMatrix(tcmd.levers[pi]);
        Eigen::RowVector3d dotCrossDiff(tcmd.axis.transpose()*crossMat);
        double squareDiff = 2.*tcmd.axis.dot(lever)*tcmd.scale*scale_;


        const auto& qDiffX = tcmd.jac.vectorJacobian(rd.pgdata->mb(), rd.pgdata->mbc(),
                                                     Eigen::Vector3d::UnitX()).block(3, 0, 3, tcmd.jac.dof());
        tcmd.jacMatTmp.row(0) = forceData.forces[pi].force().transpose()*qDiffX;

        const auto& qDiffY = tcmd.jac.vectorJacobian(rd.pgdata->mb(), rd.pgdata->mbc(),
                                                     Eigen::Vector3d::UnitY()).block(3, 0, 3, tcmd.jac.dof());
        tcmd.jacMatTmp.row(1) = forceData.forces[pi].force().transpose()*qDiffY;

        const auto& qDiffZ = tcmd.jac.vectorJacobian(rd.pgdata->mb(), rd.pgdata->mbc(),
                                                     Eigen::Vector3d::UnitZ()).block(3, 0, 3, tcmd.jac.dof());
        tcmd.jacMatTmp.row(2) = forceData.forces[pi].force().transpose()*qDiffZ;

        tcmd.jacMat.noalias() = (squareDiff*dotCrossDiff)*tcmd.jacMatTmp;

        updateFullJacobianSparse(rd.pgdata->mb(), tcmd.jac, tcmd.jacMat, tcmd.jacMatFull,
                                 {0, rd.pgdata->qParamsBegin()});
        gradient += tcmd.jacMatFull.transpose();

        Eigen::RowVector3d fDiff;
        fDiff = squareDiff*dotCrossDiff*X_0_b.rotation();

        gradient.coeffRef(forceIndex + 0) += fDiff.x();
        gradient.coeffRef(forceIndex + 1) += fDiff.y();
        gradient.coeffRef(forceIndex + 2) += fDiff.z();
        forceIndex += 3;
      }
    }


    for(NormalForceTargetData& nft: rd.normalForceTargets)
    {
      const auto& forceData = rd.pgdata->forceDatas()[nft.forcePos];
      const sva::PTransformd& X_0_b = rd.pgdata->mbc().bodyPosW[forceData.bodyIndex];

      double nForceSum = 0.;
      for(std::size_t pi = 0; pi < forceData.points.size(); ++pi)
      {
        // force in the point pi frame
        sva::PTransformd X_0_pi = forceData.points[pi]*X_0_b;
        nForceSum += (X_0_pi.rotation()*forceData.forces[pi].force()).z();
      }
      double error = nForceSum - nft.target;

      int forceIndex = int(nft.gradientPos);
      double scale = scale_*nft.scale;
      for(std::size_t pi = 0; pi < forceData.points.size(); ++pi)
      {
        double coef = 2.*scale*error;
        sva::PTransformd X_0_pi = forceData.points[pi]*X_0_b;

        const auto& jacPointMat =
            nft.jac.vectorJacobian(rd.pgdata->mb(),
                                   rd.pgdata->mbc(),
                                   forceData.points[pi].rotation().row(2).transpose())\
            .block(3, 0, 3, nft.jac.dof());
        nft.jacMat.noalias() = (coef*forceData.forces[pi].force().transpose())*
          jacPointMat;
        updateFullJacobianSparse(rd.pgdata->mb(), nft.jac, nft.jacMat,
                                 nft.jacMatFull,
                                 {0, rd.pgdata->qParamsBegin()});
        gradient += nft.jacMatFull.transpose();

        Eigen::RowVector3d fDiff;
        fDiff = coef*X_0_pi.rotation().row(2);
        gradient.coeffRef(forceIndex + 0) += fDiff.x();
        gradient.coeffRef(forceIndex + 1) += fDiff.y();
        gradient.coeffRef(forceIndex + 2) += fDiff.z();
        forceIndex += 3;
      }
    }


    for(TangentialForceMinimizationData& tfm: rd.tanForceMin)
    {
      const auto& forceData = rd.pgdata->forceDatas()[tfm.forcePos];
      const sva::PTransformd& X_0_b = rd.pgdata->mbc().bodyPosW[forceData.bodyIndex];

      int forceIndex = int(tfm.gradientPos);
      double scale = scale_*tfm.scale;
      for(std::size_t pi = 0; pi < forceData.points.size(); ++pi)
      {
        sva::PTransformd X_0_pi = forceData.points[pi]*X_0_b;
        Eigen::Vector3d fPoint(X_0_pi.rotation()*forceData.forces[pi].force());

        double coefT = 2.*scale*fPoint.x();
        double coefB = 2.*scale*fPoint.y();

        auto fillGradient =
          [forceIndex, pi, &X_0_pi, &forceData, &rd, &gradient, &tfm](int row, double coef)
        {
          const auto& jacPointMat =
              tfm.jac.vectorJacobian(rd.pgdata->mb(),
                                     rd.pgdata->mbc(),
                                     forceData.points[pi].rotation().row(row).transpose())\
              .block(3, 0, 3, tfm.jac.dof());
          tfm.jacMat.noalias() = (coef*forceData.forces[pi].force().transpose())*
            jacPointMat;
          updateFullJacobianSparse(rd.pgdata->mb(), tfm.jac, tfm.jacMat,
                                   tfm.jacMatFull,
                                   {0, rd.pgdata->qParamsBegin()});
          gradient += tfm.jacMatFull.transpose();

          Eigen::RowVector3d fDiff;
          fDiff = coef*X_0_pi.rotation().row(row);
          gradient.coeffRef(forceIndex + 0) += fDiff.x();
          gradient.coeffRef(forceIndex + 1) += fDiff.y();
          gradient.coeffRef(forceIndex + 2) += fDiff.z();
        };
        fillGradient(0, coefT);
        fillGradient(1, coefB);

        forceIndex += 3;
      }
    }


    for(BodyPositionTargetData& bp: rd.bodyPosTargets)
    {
      Eigen::Vector3d error(rd.pgdata->mbc().bodyPosW[bp.bodyIndex].translation()
        - bp.target);

      const Eigen::MatrixXd& jacMat = bp.jac.jacobian(rd.pgdata->mb(), rd.pgdata->mbc());
      bp.jacMat.noalias() = (scale_*bp.scale*2.*error.transpose())*\
          jacMat.block(3, 0, 3, bp.jac.dof());
      updateFullJacobianSparse(rd.pgdata->mb(), bp.jac, bp.jacMat, bp.jacMatFull,
                               {0, rd.pgdata->qParamsBegin()});
      gradient += bp.jacMatFull.transpose();
    }


    for(BodyOrientationTargetData& bo: rd.bodyOriTargets)
    {
      Eigen::Vector3d error(sva::rotationError(bo.target,
        rd.pgdata->mbc().bodyPosW[bo.bodyIndex].rotation(), 1e-7));

      const Eigen::MatrixXd& jacMat = bo.jac.jacobian(rd.pgdata->mb(), rd.pgdata->mbc());
      bo.jacMat.noalias() = (scale_*bo.scale*2.*error.transpose())*\
          jacMat.block(0, 0, 3, bo.jac.dof());
      updateFullJacobianSparse(rd.pgdata->mb(), bo.jac, bo.jacMat, bo.jacMatFull,
                               {0, rd.pgdata->qParamsBegin()});
      gradient += bo.jacMatFull.transpose();
    }
  }
}


} // pg
