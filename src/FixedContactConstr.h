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
// roboptim
#include <roboptim/core/differentiable-function.hh>

// RBDyn
#include <RBDyn/Jacobian.h>

// PG
#include "PGData.h"

namespace pg
{

class FixedPositionContactConstr : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  FixedPositionContactConstr(PGData* pgdata, int bodyId,
      const Eigen::Vector3d& target,
      const sva::PTransformd& surfaceFrame)
    : roboptim::DifferentiableFunction(pgdata->pbSize(), 3, "FixedPositionContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , target_(target)
    , surfaceFrame_(surfaceFrame)
    , jac_(pgdata->multibody(), bodyId, surfaceFrame.translation())
  {}
  ~FixedPositionContactConstr() throw()
  { }


  void impl_compute(result_t& res, const argument_t& x) const throw()
  {
   // std::cout << x.transpose() << std::endl;
    pgdata_->x(x);
    sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
    res = pos.translation() - target_;
    std::cout << target_.transpose() << std::endl;
    std::cout << pos.translation().transpose() << std::endl;
    std::cout << res.transpose() << std::endl;
    std::cout << std::endl;
  }


  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    const Eigen::MatrixXd& jacMat = jac_.jacobian(pgdata_->multibody(), pgdata_->mbc());
    jac_.fullJacobian(pgdata_->multibody(), jacMat.block(3, 0, 3, jacMat.cols()), jac);
    std::cout << jac << std::endl;
    std::cout << std::endl;
  }


  void impl_gradient(gradient_t& gradient,
      const argument_t& x, size_type functionId) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;

  int bodyIndex_;
  Eigen::Vector3d target_;
  sva::PTransformd surfaceFrame_;
  mutable rbd::Jacobian jac_;
};




class FixedOrientationContactConstr : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  FixedOrientationContactConstr(PGData* pgdata, int bodyId,
      const Eigen::Matrix3d& target,
      const sva::PTransformd& surfaceFrame)
    : roboptim::DifferentiableFunction(pgdata->pbSize(), 3, "FixedOrientationContact")
    , pgdata_(pgdata)
    , bodyIndex_(pgdata->multibody().bodyIndexById(bodyId))
    , target_(target)
    , surfaceFrame_(surfaceFrame)
    , jac_(pgdata->multibody(), bodyId)
  {}
  ~FixedOrientationContactConstr() throw()
  { }


  void impl_compute(result_t& res, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
    res(0) = pos.rotation().row(0).dot(target_.row(0));
    res(1) = pos.rotation().row(1).dot(target_.row(1));
    res(2) = pos.rotation().row(2).dot(target_.row(2));
  }

  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    sva::PTransformd pos = surfaceFrame_*pgdata_->mbc().bodyPosW[bodyIndex_];
    {
      const Eigen::MatrixXd& mat = jac_.vectorBodyJacobian(pgdata_->multibody(), pgdata_->mbc(), pos.rotation().row(0).transpose());
      //Eigen::MatrixXd tmp = ((mat.array().colwise()*target_.row(0).transpose().array()).array().rowwise().sum()).matrix();
      Eigen::MatrixXd tmp = mat.transpose()*target_.row(0).transpose();
      Eigen::MatrixXd tmp2(1, pgdata_->multibody().nrDof());
      jac_.fullJacobian(pgdata_->multibody(), tmp, tmp2);
      jac.row(0) = tmp2;
    }
    {
      const Eigen::MatrixXd& mat = jac_.vectorBodyJacobian(pgdata_->multibody(), pgdata_->mbc(), pos.rotation().row(1)).transpose();
      //Eigen::MatrixXd tmp = ((mat.array().colwise()*target_.row(1).transpose().array()).array().rowwise().sum()).matrix();
      Eigen::MatrixXd tmp = mat.transpose()*target_.row(1).transpose();
      Eigen::MatrixXd tmp2(1, pgdata_->multibody().nrDof());
      jac_.fullJacobian(pgdata_->multibody(), tmp, tmp2);
      jac.row(1) = tmp2;
    }
    {
      const Eigen::MatrixXd& mat = jac_.vectorBodyJacobian(pgdata_->multibody(), pgdata_->mbc(), pos.rotation().row(2).transpose());
      //Eigen::MatrixXd tmp = ((mat.array().colwise()*target_.row(2).transpose().array()).array().rowwise().sum()).matrix();
      Eigen::MatrixXd tmp = mat.transpose()*target_.row(2).transpose();
      Eigen::MatrixXd tmp2(1, pgdata_->multibody().nrDof());
      jac_.fullJacobian(pgdata_->multibody(), tmp, tmp2);
      jac.row(2) = tmp2;
    }
  }

  void impl_gradient(gradient_t& gradient,
      const argument_t& x, size_type functionId) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;

  int bodyIndex_;
  sva::Matrix3d target_;
  sva::PTransformd surfaceFrame_;
  mutable rbd::Jacobian jac_;
};

} // namespace pg
