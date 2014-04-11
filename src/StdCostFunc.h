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
// boost
#include <boost/math/constants/constants.hpp>

// roboptim
#include <roboptim/core/differentiable-function.hh>

// PG
#include "ConfigStruct.h"
#include "PGData.h"

namespace pg
{

class StdCostFunc : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  StdCostFunc(PGData* pgdata, std::vector<std::vector<double>> q,
              double postureScale, double torqueScale, double forceScale, 
              double ellipseScale,
              const std::vector<BodyPositionTarget>& bodyPosTargets,
              const std::vector<BodyOrientationTarget>& bodyOriTargets,
              const std::vector<ForceContact>& forceContacts,
              const std::vector<ForceContactMinimization>& forceContactsMin)
    : roboptim::DifferentiableFunction(pgdata->pbSize(), 1, "StdCostFunc")
    , pgdata_(pgdata)
    , tq_(std::move(q))
    , postureScale_(postureScale)
    , torqueScale_(torqueScale)
    , forceScale_(forceScale)
    , ellipseScale_(ellipseScale)
    , bodyPosTargets_(bodyPosTargets.size())
    , bodyOriTargets_(bodyOriTargets.size())
    , forceContactsMin_()
  {
    for(std::size_t i = 0; i < bodyPosTargets_.size(); ++i)
    {
      rbd::Jacobian jac(pgdata->mb(), bodyPosTargets[i].bodyId);
      Eigen::MatrixXd jacMat(1, jac.dof());
      Eigen::MatrixXd jacMatFull(1, pgdata->mb().nrDof());
      bodyPosTargets_[i] = {pgdata->multibody().bodyIndexById(bodyPosTargets[i].bodyId),
                            bodyPosTargets[i].target,
                            bodyPosTargets[i].scale,
                            jac, jacMat, jacMatFull};
    }

    for(std::size_t i = 0; i < bodyOriTargets_.size(); ++i)
    {
      rbd::Jacobian jac(pgdata->mb(), bodyPosTargets[i].bodyId);
      Eigen::MatrixXd jacMat(1, jac.dof());
      Eigen::MatrixXd jacMatFull(1, pgdata->mb().nrDof());
      bodyOriTargets_[i] = {pgdata->multibody().bodyIndexById(bodyOriTargets[i].bodyId),
                            bodyOriTargets[i].target,
                            bodyOriTargets[i].scale,
                            jac, jacMat, jacMatFull};
    }

    for(std::size_t i = 0; i < forceContactsMin.size(); ++i)
    {
      std::size_t gradientPos = pgdata->forceParamsBegin();
      for(std::size_t j = 0; j < forceContacts.size(); ++j)
      {
        // we don't break scine it could be many contact on the same body
        if(forceContactsMin[i].bodyId == forceContacts[j].bodyId)
        {
          forceContactsMin_.push_back({j, gradientPos,forceContactsMin[i].scale});
        }
        gradientPos += forceContacts[j].points.size()*3;
      }
    }
  }


  void impl_compute(result_t& res, const argument_t& x) const throw()
  {
    pgdata_->x(x);

    // compute posture task
    double posture = 0.;
    if(postureScale_ > 0.)
    {
      const std::vector<std::vector<double>>& q = pgdata_->mbc().q;
      for(int i = 0; i < pgdata_->multibody().nrJoints(); ++i)
      {
        if(pgdata_->multibody().joint(i).params() == 1)
        {
          posture += std::pow(q[i][0] - tq_[i][0], 2);
        }
      }
    }

    // compute torque task
    /*
    scalar_t torque = scalar_t(0., Eigen::VectorXd::Zero(this->inputSize()));
    if(torqueScale_ > 0.)
    {
      const auto& torqueVec = pgdata_->id().torque();
      for(int i = 0; i < pgdata_->multibody().nrJoints(); ++i)
      {
        for(int j = 0; j < pgdata_->multibody().joint(i).dof(); ++j)
        {
          torque += std::pow(torqueVec[i](j), 2);
        }
      }
    }
    */

    // compute force task
    double force = 0.;
    Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
    if(forceScale_ > 0.)
    {
      for(const auto& fd: pgdata_->forceDatas())
      {
        Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
        for(const sva::ForceVecd& fv: fd.forces)
        {
          forceTmp += fv.force();
        }
        force += forceTmp.squaredNorm()*forceScale_;
      }
    }

    double pos = 0.;
    for(const BodyPositionTargetData& bp: bodyPosTargets_)
    {
      pos += (pgdata_->mbc().bodyPosW[bp.bodyIndex].translation() - bp.target).squaredNorm()*
          bp.scale;
    }

    double ori = 0.;
    for(const BodyOrientationTargetData& bo: bodyOriTargets_)
    {
      ori += sva::rotationError(bo.target,
        pgdata_->mbc().bodyPosW[bo.bodyIndex].rotation(), 1e-7).squaredNorm()*
          bo.scale;
    }

    double forceMin = 0.;
    for(const ForceContactMinimizationData& fcmd: forceContactsMin_)
    {
      const auto& forceData = pgdata_->forceDatas()[fcmd.forcePos];
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
      for(std::size_t i = 0; i < pgdata_->ellipseDatas().size(); ++i)
      {
        ellipses += -boost::math::constants::pi<double>()*pgdata_->ellipseDatas()[i].r1 * pgdata_->ellipseDatas()[i].r2; 
        std::cout << "ellipses cost: " << ellipses << std::endl;
      } 
    }
    */

    res(0) = posture*postureScale_ + force +
        pos + ori + forceMin;
  }

  void impl_gradient(gradient_t& gradient,
      const argument_t& x, size_type /* functionId */) const throw()
  {
    pgdata_->x(x);
    gradient.setZero();

    if(postureScale_ > 0.)
    {
      int index = 0;
      const std::vector<std::vector<double>>& q = pgdata_->mbc().q;
      double coef = 2.*postureScale_;
      for(int i = 0; i < pgdata_->multibody().nrJoints(); ++i)
      {
        if(pgdata_->multibody().joint(i).params() == 1)
        {
          gradient(index) += coef*(q[i][0] - tq_[i][0]);
        }
        index += pgdata_->mb().joint(i).dof();
      }
    }


    if(forceScale_ > 0.)
    {
      int index = pgdata_->forceParamsBegin();
      for(const auto& fd: pgdata_->forceDatas())
      {
        Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
        for(const sva::ForceVecd& fv: fd.forces)
        {
          forceTmp += fv.force();
        }
        forceTmp *= 2.*forceScale_;

        for(std::size_t i = 0; i < fd.forces.size(); ++i)
        {
          gradient(index + 0) += forceTmp.x();
          gradient(index + 1) += forceTmp.y();
          gradient(index + 2) += forceTmp.z();
          index += 3;
        }
      }
    }


    for(const ForceContactMinimizationData& fcmd: forceContactsMin_)
    {
      const auto& forceData = pgdata_->forceDatas()[fcmd.forcePos];
      Eigen::Vector3d forceTmp(Eigen::Vector3d::Zero());
      std::size_t index = fcmd.gradientPos;

      for(const sva::ForceVecd& fv: forceData.forces)
      {
        forceTmp += fv.force();
      }
      forceTmp *= 2.*fcmd.scale;

      for(std::size_t i = 0; i < forceData.forces.size(); ++i)
      {
        gradient(index + 0) += forceTmp.x();
        gradient(index + 1) += forceTmp.y();
        gradient(index + 2) += forceTmp.z();
        index += 3;
      }
    }


    for(BodyPositionTargetData& bp: bodyPosTargets_)
    {
      Eigen::Vector3d error(pgdata_->mbc().bodyPosW[bp.bodyIndex].translation()
        - bp.target);

      const Eigen::MatrixXd& jacMat = bp.jac.jacobian(pgdata_->mb(), pgdata_->mbc());
      bp.jacMat.noalias() = (bp.scale*2.*error.transpose())*jacMat.block(3, 0, 3, bp.jac.dof());
      bp.jac.fullJacobian(pgdata_->mb(), bp.jacMat, bp.jacMatFull);
      gradient.head(bp.jacMatFull.cols()).noalias() += bp.jacMatFull.transpose();
    }


    for(BodyOrientationTargetData& bo: bodyOriTargets_)
    {
      Eigen::Vector3d error(sva::rotationError(bo.target,
        pgdata_->mbc().bodyPosW[bo.bodyIndex].rotation(), 1e-7));

      const Eigen::MatrixXd& jacMat = bo.jac.jacobian(pgdata_->mb(), pgdata_->mbc());
      bo.jacMat.noalias() = (bo.scale*2.*error.transpose())*jacMat.block(0, 0, 3, bo.jac.dof());
      bo.jac.fullJacobian(pgdata_->mb(), bo.jacMat, bo.jacMatFull);
      gradient.head(bo.jacMatFull.cols()).noalias() += bo.jacMatFull.transpose();
    }
  }

private:
  struct BodyPositionTargetData
  {
    int bodyIndex;
    Eigen::Vector3d target;
    double scale;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat, jacMatFull;
  };

  struct BodyOrientationTargetData
  {
    int bodyIndex;
    Eigen::Matrix3d target;
    double scale;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat, jacMatFull;
  };

  struct ForceContactMinimizationData
  {
    std::size_t forcePos;
    std::size_t gradientPos;
    double scale;
  };

private:
  PGData* pgdata_;
  std::vector<std::vector<double>> tq_;
  double postureScale_;
  double torqueScale_;
  double forceScale_;
  double ellipseScale_;
  mutable std::vector<BodyPositionTargetData> bodyPosTargets_;
  mutable std::vector<BodyOrientationTargetData> bodyOriTargets_;
  std::vector<ForceContactMinimizationData> forceContactsMin_;
};

} // namespace pg
