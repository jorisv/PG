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
#include "AutoDiffFunction.h"
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


template<typename T>
T computeDist(SCD::CD_Pair* pair,
              sva::PTransform<T> obj1, sva::PTransform<T> obj2,
              sva::PTransformd obj1d, sva::PTransformd obj2d)
{
  using namespace Eigen;

  SCD::Point3 pb1Tmp, pb2Tmp;
  double dist = pair->getClosestPoints(pb1Tmp, pb2Tmp);
  double sign = dist >= 0. ? 1.: -1.;

  sva::PTransformd pb1World(Vector3d(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]));
  sva::PTransformd pb2World(Vector3d(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]));

  Vector3d pb1Local = (pb1World*obj1d.inv()).translation();
  Vector3d pb2Local = (pb2World*obj2d.inv()).translation();

  Vector3<T> pb1((sva::PTransform<T>(Vector3<T>(pb1Local.cast<T>()))*obj1).translation());
  Vector3<T> pb2((sva::PTransform<T>(Vector3<T>(pb2Local.cast<T>()))*obj2).translation());

  return sign*(pb2 - pb1).squaredNorm();
}


template<typename Type>
class EnvCollisionConstr : public AutoDiffFunction<Type, Eigen::Dynamic>
{
public:
  typedef AutoDiffFunction<Type, Eigen::Dynamic> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  EnvCollisionConstr(PGData<Type>* pgdata, const std::vector<EnvCollision>& cols)
    : parent_t(pgdata->pbSize(), int(cols.size()), "EnvCollision")
    , pgdata_(pgdata)
  {
    cols_.reserve(cols.size());
    for(const EnvCollision& sc: cols)
    {
      cols_.push_back({pgdata_->multibody().bodyIndexById(sc.bodyId),
                       sc.bodyT, new SCD::CD_Pair(sc.bodyHull, sc.envHull)});
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


  void impl_compute(result_ad_t& res, const argument_t& x) const
  {
    pgdata_->x(x);
    const FK<scalar_t>& fk = pgdata_->fk();
    int i = 0;
    for(const CollisionData& cd: cols_)
    {
      const sva::PTransform<scalar_t>& objPos = fk.bodyPosW()[cd.bodyIndex];
      sva::PTransformd objPosd(toValue(objPos.rotation()), toValue(objPos.translation()));

      cd.pair->operator[](0)->setTransformation(toSCD(cd.bodyT*objPosd));

      res(i) = computeDist(cd.pair,
                           objPos, sva::PTransform<scalar_t>::Identity(),
                           objPosd, sva::PTransformd::Identity());
      ++i;
    }
  }

private:
  struct CollisionData
  {
    int bodyIndex;
    sva::PTransformd bodyT;
    SCD::CD_Pair* pair;
  };

private:
  PGData<Type>* pgdata_;
  std::vector<CollisionData> cols_;
};

} // namespace pg

