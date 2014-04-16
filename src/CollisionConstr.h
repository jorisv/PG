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

// SCD
#include <SCD/Matrix/SCD_Types.h>


// forward declaration
namespace SCD
{
class CD_Pair;
} // SCD

namespace pg
{
class PGData;
class EnvCollision;
class SelfCollision;


SCD::Matrix4x4 toSCD(const sva::PTransformd& t);


class EnvCollisionConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  EnvCollisionConstr(PGData* pgdata, const std::vector<EnvCollision>& cols);
  ~EnvCollisionConstr() throw();


  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  struct CollisionData
  {
    int bodyIndex;
    sva::PTransformd bodyT;
    SCD::CD_Pair* pair;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
  };

private:
  PGData* pgdata_;
  int nrNonZero_;
  mutable std::vector<CollisionData> cols_;
};



class SelfCollisionConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  SelfCollisionConstr(PGData* pgdata, const std::vector<SelfCollision>& cols);
  ~SelfCollisionConstr() throw();

  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  struct CollisionData
  {
    int body1Index;
    sva::PTransformd body1T;
    rbd::Jacobian jac1;
    Eigen::MatrixXd jac1Mat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jac1MatFull;
    int body2Index;
    sva::PTransformd body2T;
    rbd::Jacobian jac2;
    Eigen::MatrixXd jac2Mat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jac2MatFull;
    SCD::CD_Pair* pair;
  };

private:
  PGData* pgdata_;
  int nrNonZero_;
  mutable std::vector<CollisionData> cols_;
};

} // namespace pg

