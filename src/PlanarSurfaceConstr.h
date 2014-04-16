
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
#include <roboptim/core.hh>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/Jacobian.h>


namespace pg
{
class PGData;

class PlanarPositionContactConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarPositionContactConstr(PGData* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const sva::PTransformd& surfaceFrame);
  ~PlanarPositionContactConstr() throw();


  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;

  int bodyIndex_;
  sva::PTransformd targetFrame_;
  sva::PTransformd surfaceFrame_;
  mutable rbd::Jacobian jac_;
  mutable Eigen::MatrixXd dotCache_;
};



class PlanarOrientationContactConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarOrientationContactConstr(PGData* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const sva::PTransformd& surfaceFrame,
      int axis);
  ~PlanarOrientationContactConstr() throw();


  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;

  int bodyIndex_;
  sva::PTransformd targetFrame_;
  sva::PTransformd surfaceFrame_;
  int axis_;
  mutable rbd::Jacobian jac_;
  mutable Eigen::MatrixXd dotCache_;
};



class PlanarInclusionConstr : public roboptim::DifferentiableSparseFunction

{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  PlanarInclusionConstr(PGData* pgdata, int bodyId,
      const sva::PTransformd& targetFrame,
      const std::vector<Eigen::Vector2d>& targetPoints,
      const sva::PTransformd& surfaceFrame,
      const std::vector<Eigen::Vector2d>& surfacePoints);
  ~PlanarInclusionConstr() throw();


  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;

  int bodyIndex_;
  sva::PTransformd targetFrame_;
  std::vector<Eigen::Vector2d> targetPoints_;
  std::vector<Eigen::Vector2d> targetVecNorm_;
  std::vector<sva::PTransformd> surfacePoints_; ///< Surface points in body coord.
  mutable rbd::Jacobian jac_;
  mutable Eigen::MatrixXd transJac_;
  mutable Eigen::MatrixXd tJac_;
  mutable Eigen::MatrixXd bJac_;
  mutable Eigen::MatrixXd sumJac_;
  mutable Eigen::MatrixXd fullJac_;
};

} // namespace pg
