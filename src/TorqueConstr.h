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
// std
#include <set>

// Eigen
#include <unsupported/Eigen/Polynomials>

// PG
#include "AutoDiffFunction.h"
#include "PGData.h"

namespace pg
{

template<typename Type>
class TorqueConstr : public AutoDiffFunction<Type, Eigen::Dynamic>
{
public:
  typedef AutoDiffFunction<Type, Eigen::Dynamic> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;
  typedef typename parent_t::intervals_t intervals_t;

public:
  TorqueConstr(PGData<Type>* pgdata, const std::vector<SpringJoint>& springs,
               const Eigen::VectorXd& tl, const Eigen::VectorXd& tu)
    : parent_t(pgdata, pgdata->pbSize(),
               (pgdata->multibody().nrDof() - pgdata->multibody().joint(0).dof()),
               "Torque")
    , pgdata_(pgdata)
  {
    // store spring joint index to avoid to add a torque constraint on them
    std::set<int> springJoints;
    for(const SpringJoint& sj: springs)
    {
      int index = pgdata->multibody().jointIndexById(sj.jointId);
      assert(pgdata->multibody().joint(index).dof() == 1);
      springJ_.push_back({index, sj.K, sj.O});
      springJoints.insert(index);
      bounds_.push_back({0., 0.});
    }

    // don't take the joint 0 (free flyier or fixed joint)
    for(int i = 1; i < pgdata->multibody().nrJoints(); ++i)
    {
      // only add it if it's not a spring
      if(springJoints.find(i) == springJoints.end())
      {
        actJ_.push_back({i});
        int posInParam = pgdata->multibody().jointPosInDof(i) - pgdata->multibody().joint(0).dof();
        for(int dof = 0; dof < pgdata->multibody().joint(i).dof(); ++dof)
        {
          bounds_.push_back({tl[posInParam + dof], tu[posInParam + dof]});
        }
      }
    }
  }
  ~TorqueConstr() throw()
  { }


  intervals_t bounds() const
  {
    return bounds_;
  }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
    const ID<scalar_t>& id = pgdata_->id();

    int pos = 0;
    for(const SpringJointData& sj: springJ_)
    {
      const auto& t = id.torque()[sj.index];
      scalar_t spring = sj.K*pgdata_->q()[sj.index][0] + sj.O;
      res(pos) = t(0) - spring;
      ++pos;
    }

    for(const ActuatedJointData& aj: actJ_)
    {
      const auto& t = id.torque()[aj.index];
      res.segment(pos, t.size()) = t;
      pos += int(t.size());
    }
  }

private:
  struct ActuatedJointData
  {
    int index;
  };

  struct SpringJointData
  {
    int index;
    double K;
    double O;
  };

private:
  PGData<Type>* pgdata_;
  std::vector<ActuatedJointData> actJ_;
  std::vector<SpringJointData> springJ_;
  intervals_t bounds_;
};


template<typename Type>
class TorquePolyBoundsConstr : public AutoDiffFunction<Type, Eigen::Dynamic>
{
public:
  typedef AutoDiffFunction<Type, Eigen::Dynamic> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  TorquePolyBoundsConstr(PGData<Type>* pgdata,
                         std::vector<std::vector<Eigen::VectorXd>> tl,
                         std::vector<std::vector<Eigen::VectorXd>> tu)
    : parent_t(pgdata, pgdata->pbSize(),
               (pgdata->multibody().nrDof() - pgdata->multibody().joint(0).dof())*2,
               "TorquePolyBounds")
    , pgdata_(pgdata)
    , tl_(std::move(tl))
    , tu_(std::move(tu))
  {}
  ~TorquePolyBoundsConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
    const ID<scalar_t>& id = pgdata_->id();
    const std::vector<std::vector<scalar_t>>&q = pgdata_->q();
    int vecPos = 0;
    int dof = (pgdata_->multibody().nrDof() - pgdata_->multibody().joint(0).dof());
    for(std::size_t i = 1; i < q.size(); ++i)
    {
      for(std::size_t j = 0; j < q[i].size(); ++j)
      {
        res(vecPos) = id.torque()[i][j] - Eigen::poly_eval(tl_[i][j], q[i][j]);
        res(vecPos + dof) = id.torque()[i][j] - Eigen::poly_eval(tu_[i][j], q[i][j]);
        ++vecPos;
      }
    }
  }

private:
  PGData<Type>* pgdata_;
  std::vector<std::vector<Eigen::VectorXd>> tl_;
  std::vector<std::vector<Eigen::VectorXd>> tu_;
};

} // namespace pg
