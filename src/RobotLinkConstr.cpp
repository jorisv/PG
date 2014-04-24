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
#include "RobotLinkConstr.h"

// include
// PG
#include "PGData.h"
#include "FillSparse.h"


namespace pg
{


/*
 *                 RobotLinkConstr
 */


RobotLinkConstr::RobotLinkConstr(PGData* pgdata1, PGData* pgdata2,
      const std::vector<int>& bodiesId)
  : roboptim::DifferentiableSparseFunction(pgdata1->pbSize(), int(6*bodiesId.size()),
                                           "RobotLink")
  , pgdata1_(pgdata1)
  , pgdata2_(pgdata2)
  , nrNonZero_(0)
  , links_()
{
  for(int bodyId: bodiesId)
  {
    rbd::Jacobian jac1(pgdata1_->mb(), bodyId);
    rbd::Jacobian jac2(pgdata2_->mb(), bodyId);
    Eigen::MatrixXd jacMat1(6, jac1.dof());
    Eigen::MatrixXd jacMat2(6, jac2.dof());

    links_.push_back({jac1, jac2, jacMat1, jacMat2});

    nrNonZero_ += jac1.dof()*6 + jac2.dof()*6;
  }
}


RobotLinkConstr::~RobotLinkConstr() throw()
{ }


void RobotLinkConstr::impl_compute(result_t& res, const argument_t& x) const throw()
{
  pgdata1_->x(x);
  pgdata2_->x(x);

  for(std::size_t i = 0; i < links_.size(); ++i)
  {
    LinkData& link = links_[i];
    int index1 = link.jac1.jointsPath().back();
    int index2 = link.jac2.jointsPath().back();

    const sva::PTransformd body1 = pgdata1_->mbc().bodyPosW[index1];
    const sva::PTransformd body2 = pgdata2_->mbc().bodyPosW[index2];

    Eigen::Vector3d posErr = body1.translation() - body2.translation();
    /// @todo try to use rotation error
    Eigen::Vector3d rotErr;
    rotErr(0) = body1.rotation().row(0).dot(body2.rotation().row(0));
    rotErr(1) = body1.rotation().row(1).dot(body2.rotation().row(1));
    rotErr(2) = body1.rotation().row(2).dot(body2.rotation().row(2));

    res.segment<3>(i*6) = rotErr;
    res.segment<3>(i*6 + 3) = posErr;
  }
}


void RobotLinkConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
{
  jac.reserve(nrNonZero_);
  jac.setZero();

  pgdata1_->x(x);
  pgdata2_->x(x);

  for(std::size_t i = 0; i < links_.size(); ++i)
  {
    LinkData& link = links_[i];
    int index1 = link.jac1.jointsPath().back();
    int index2 = link.jac2.jointsPath().back();

    const sva::PTransformd body1 = pgdata1_->mbc().bodyPosW[index1];
    const sva::PTransformd body2 = pgdata2_->mbc().bodyPosW[index2];

    link.jacMat1.block(3, 0, 3, link.jac1.dof()).noalias() =
      link.jac1.jacobian(pgdata1_->mb(), pgdata1_->mbc()).block(3, 0, 3, link.jac1.dof());
    link.jacMat2.block(3, 0, 3, link.jac2.dof()).noalias() =
      -link.jac2.jacobian(pgdata2_->mb(), pgdata2_->mbc()).block(3, 0, 3, link.jac2.dof());

    for(int j = 0; j < 3; ++j)
    {
      Eigen::Vector3d axis(Eigen::Vector3d::Zero());
      axis(j) = 1.;

      const Eigen::MatrixXd& mat1 =
          link.jac1.vectorJacobian(pgdata1_->multibody(), pgdata1_->mbc(), axis);
      link.jacMat1.row(j).noalias() = body2.rotation().row(j)*mat1.block(3, 0, 3, mat1.cols());
      const Eigen::MatrixXd& mat2 =
          link.jac2.vectorJacobian(pgdata2_->multibody(), pgdata2_->mbc(), axis);
      link.jacMat2.row(j).noalias() = body1.rotation().row(j)*mat2.block(3, 0, 3, mat2.cols());
    }

    incrementFullJacobianSparse(pgdata1_->mb(), link.jac1, link.jacMat1, jac,
                             {i*6, pgdata1_->qParamsBegin()});
    incrementFullJacobianSparse(pgdata2_->mb(), link.jac2, link.jacMat2, jac,
                             {i*6, pgdata2_->qParamsBegin()});
  }
}


} // pg

