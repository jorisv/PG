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
// SCD
#include <SCD/CD/CD_Pair.h>

// PG
#include "ConfigStruct.h"
#include "PGData.h"

namespace pg
{


SCD::Matrix4x4 toSCD(const sva::PTransformd& t)
{
  SCD::Matrix4x4 m;
  const Eigen::Matrix3d& rot = t.rotation();
  const Eigen::Vector3d& tran = t.translation();

  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      m(i,j) = rot(j,i);
    }
  }

  m(0,3) = tran(0);
  m(1,3) = tran(1);
  m(2,3) = tran(2);

  return m;
}


// return the pair signed squared distance
double distance(SCD::CD_Pair* pair)
{
  return pair->getDistance();
}


// return the pair signed squared distance and the closest point position in world frame
std::tuple<double, Eigen::Vector3d, Eigen::Vector3d>
closestPoints(SCD::CD_Pair* pair)
{
  using namespace Eigen;

  SCD::Point3 pb1Tmp, pb2Tmp;
  double dist = pair->getClosestPoints(pb1Tmp, pb2Tmp);

  Vector3d T_0_p1(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]);
  Vector3d T_0_p2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

  return std::make_tuple(dist, T_0_p1, T_0_p2);
}



class EnvCollisionConstr : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  EnvCollisionConstr(PGData* pgdata, const std::vector<EnvCollision>& cols)
    : roboptim::DifferentiableFunction(pgdata->pbSize(), int(cols.size()), "EnvCollision")
    , pgdata_(pgdata)
  {
    cols_.reserve(cols.size());
    for(const EnvCollision& sc: cols)
    {
      rbd::Jacobian jac(pgdata_->mb(), sc.bodyId);
      Eigen::MatrixXd jacMat(1, jac.dof());
      Eigen::MatrixXd jacMatFull(1, pgdata_->mb().nrDof());
      cols_.push_back({pgdata_->multibody().bodyIndexById(sc.bodyId),
                       sc.bodyT, new SCD::CD_Pair(sc.bodyHull, sc.envHull),
                       jac, jacMat, jacMatFull});
    }
  }
  ~EnvCollisionConstr() throw()
  {
    /// @todo try to use unique_ptr instead
    for(auto& cd: cols_)
    {
      delete cd.pair;
    }
  }


  void impl_compute(result_t& res, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    int i = 0;
    for(const CollisionData& cd: cols_)
    {
      sva::PTransformd X_0_b(pgdata_->mbc().bodyPosW[cd.bodyIndex]);

      cd.pair->operator[](0)->setTransformation(toSCD(cd.bodyT*X_0_b));

      res(i) = distance(cd.pair);
      ++i;
    }
  }

  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw()
  {
    pgdata_->x(x);
    jac.setZero();

    int i = 0;
    for(CollisionData& cd: cols_)
    {
      sva::PTransformd X_0_b(pgdata_->mbc().bodyPosW[cd.bodyIndex]);

      cd.pair->operator[](0)->setTransformation(toSCD(cd.bodyT*X_0_b));

      double dist;
      Eigen::Vector3d T_0_p, pEnv;
      std::tie(dist, T_0_p, pEnv) = closestPoints(cd.pair);

      Eigen::Vector3d dist3d(T_0_p - pEnv);
      Eigen::Vector3d T_b_p(X_0_b.rotation()*(T_0_p - X_0_b.translation()));
      cd.jac.point(T_b_p);

      double coef = std::copysign(2., dist);
      const Eigen::MatrixXd& jacMat = cd.jac.jacobian(pgdata_->mb(), pgdata_->mbc());
      cd.jacMat.noalias() = coef*dist3d.transpose()*jacMat.block(3, 0, 3, cd.jac.dof());
      cd.jac.fullJacobian(pgdata_->mb(), cd.jacMat, cd.jacMatFull);
      jac.row(i) = cd.jacMatFull;
      ++i;
    }
  }

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
    Eigen::MatrixXd jacMat, jacMatFull;
  };

private:
  PGData* pgdata_;
  mutable std::vector<CollisionData> cols_;
};



/*
template<typename Type>
class SelfCollisionConstr : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  SelfCollisionConstr(PGData* pgdata, const std::vector<SelfCollision>& cols)
    : parent_t(pgdata->pbSize(), int(cols.size()), "EnvCollision")
    , pgdata_(pgdata)
  {
    cols_.reserve(cols.size());
    for(const SelfCollision& sc: cols)
    {
      cols_.push_back({pgdata_->multibody().bodyIndexById(sc.body1Id),
                       sc.body1T,
                       pgdata_->multibody().bodyIndexById(sc.body2Id),
                       sc.body2T,
                       new SCD::CD_Pair(sc.body1Hull, sc.body2Hull)});
    }
  }
  ~SelfCollisionConstr() throw()
  {
    /// @todo try to use unique_ptr instead
    for(auto& cd: cols_)
    {
      delete cd.pair;
    }
  }


  void impl_compute(result_ad_t& res, const argument_t& x) const
  {
    const FK<scalar_t>& fk = pgdata_->fk();
    int i = 0;
    for(const CollisionData& cd: cols_)
    {
      const sva::PTransform<scalar_t>& obj1Pos = fk.bodyPosW()[cd.body1Index];
      const sva::PTransform<scalar_t>& obj2Pos = fk.bodyPosW()[cd.body2Index];
      sva::PTransformd obj1Posd(toValue(obj1Pos.rotation()), toValue(obj1Pos.translation()));
      sva::PTransformd obj2Posd(toValue(obj2Pos.rotation()), toValue(obj2Pos.translation()));

      cd.pair->operator[](0)->setTransformation(toSCD(cd.body1T*obj1Posd));
      cd.pair->operator[](1)->setTransformation(toSCD(cd.body2T*obj2Posd));

      res(i) = computeDist(cd.pair,
                           obj1Pos, obj2Pos,
                           obj1Posd, obj2Posd);
      ++i;
    }
  }

private:
  struct CollisionData
  {
    int body1Index;
    sva::PTransformd body1T;
    int body2Index;
    sva::PTransformd body2T;
    SCD::CD_Pair* pair;
  };

private:
  PGData<Type>* pgdata_;
  std::vector<CollisionData> cols_;
};
*/

} // namespace pg

