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
#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-factory.hh>
#include <roboptim/core/plugin/ipopt.hh>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// PG
#include "PGData.h"
#include "StdCostFunc.h"
#include "FixedContactConstr.h"
#include "StaticStabilityConstr.h"

namespace pg
{

struct FixedPositionContact
{
  int bodyId;
  Eigen::Vector3d target;
  sva::PTransformd surfaceFrame;
};


struct FixedOrientationContact
{
  int bodyId;
  Eigen::Matrix3d target;
  sva::PTransformd surfaceFrame;
};


struct ForceContact
{
  int bodyId;
  std::vector<sva::PTransformd> points;
  /// @todo add friction stuff
};


template<typename Type>
class PostureGenerator
{
public:
  typedef typename Type::scalar_t scalar_t;
  typedef roboptim::IpoptSolver::solver_t solver_t;

public:
  PostureGenerator(const rbd::MultiBody& mb);

  void fixedPositionContacts(std::vector<FixedPositionContact> contacts);
  void fixedOrientationContacts(std::vector<FixedOrientationContact> contacts);
  void forceContacts(std::vector<ForceContact> contacts);
  void qBounds(const std::vector<std::vector<double>>& lq,
               const std::vector<std::vector<double>>& uq);

  void param(const std::string& name, const std::string& value);
  void param(const std::string& name, double value);
  void param(const std::string& name, int value);

  bool run(const std::vector<std::vector<double>>& q);

  std::vector<std::vector<double>> q() const;

private:
  PGData<Type> pgdata_;

  solver_t::parameters_t params_;

  std::vector<FixedPositionContact> fixedPosContacts_;
  std::vector<FixedOrientationContact> fixedOriContacts_;
  std::vector<ForceContact> forceContacts_;
  Eigen::VectorXd ql_, qu_;

  Eigen::VectorXd x_;
};


// inline


/// Visitor to take the x component of the result
struct ResultVisitor : public boost::static_visitor<>
{
  void operator()(const roboptim::Result& res)
  {
    x = res.x;
  }

  void operator()(const roboptim::ResultWithWarnings& res)
  {
    x = res.x;
  }

  template <typename T>
  void operator()(const T& /* res */)
  {
    throw std::runtime_error("This result type is not handled !");
  }

  Eigen::VectorXd x;
};


template<typename Type>
PostureGenerator<Type>::PostureGenerator(const rbd::MultiBody& mb)
  : pgdata_(mb)
  , ql_(mb.nrParams())
  , qu_(mb.nrParams())
{
  ql_.setConstant(-std::numeric_limits<double>::infinity());
  qu_.setConstant(std::numeric_limits<double>::infinity());
}


template<typename Type>
void PostureGenerator<Type>::fixedPositionContacts(std::vector<FixedPositionContact> contacts)
{
  fixedPosContacts_ = std::move(contacts);
}


template<typename Type>
void PostureGenerator<Type>::fixedOrientationContacts(std::vector<FixedOrientationContact> contacts)
{
  fixedOriContacts_ = std::move(contacts);
}


template<typename Type>
void PostureGenerator<Type>::forceContacts(std::vector<ForceContact> contacts)
{
  typedef PGData<Type> pgdata_t;
  typedef typename pgdata_t::ForceData forcedata_t;

  forceContacts_ = std::move(contacts);
  std::vector<forcedata_t> forceDatas;
  forceDatas.reserve(forceContacts_.size());
  for(const ForceContact& fc: forceContacts_)
  {
    std::vector<sva::PTransform<scalar_t>> points(fc.points.size());
    std::vector<sva::ForceVec<scalar_t>> forces(fc.points.size());
    for(std::size_t i = 0; i < fc.points.size(); ++i)
    {
      points[i] = fc.points[i].cast<scalar_t>();
      forces[i] = sva::ForceVec<scalar_t>(Eigen::Vector6<scalar_t>::Zero());
    }
    forceDatas.push_back({pgdata_.multibody().bodyIndexById(fc.bodyId), points, forces});
  }
  pgdata_.forces(forceDatas);
}


template<typename Type>
void PostureGenerator<Type>::qBounds(const std::vector<std::vector<double>>& ql,
    const std::vector<std::vector<double>>& qu)
{
  ql_ = rbd::paramToVector(pgdata_.multibody(), ql);
  qu_ = rbd::paramToVector(pgdata_.multibody(), qu);
}


template<typename Type>
void PostureGenerator<Type>::param(const std::string& name, const std::string& value)
{
  params_[name].value = value;
}


template<typename Type>
void PostureGenerator<Type>::param(const std::string& name, double value)
{
  params_[name].value = value;
}


template<typename Type>
void PostureGenerator<Type>::param(const std::string& name, int value)
{
  params_[name].value = value;
}


template<typename Type>
bool PostureGenerator<Type>::run(const std::vector<std::vector<double> >& q)
{
  pgdata_.update();

  StdCostFunc<Type> cost(&pgdata_, q);

  solver_t::problem_t problem(cost);
  problem.startingPoint() = Eigen::VectorXd::Zero(pgdata_.pbSize());
  problem.startingPoint()->head(pgdata_.multibody().nrParams()) =
      rbd::paramToVector(pgdata_.multibody(), q);

  for(int i = 0; i < pgdata_.multibody().nrParams(); ++i)
  {
    problem.argumentBounds()[i] = {ql_[i], qu_[i]};
  }

  for(const FixedPositionContact& fc: fixedPosContacts_)
  {
    boost::shared_ptr<FixedPositionContactConstr<Type>> fcc(
        new FixedPositionContactConstr<Type>(&pgdata_, fc.bodyId, fc.target, fc.surfaceFrame));
    problem.addConstraint(fcc, {{0., 0.}, {0., 0.}, {0., 0.}},
        {{1.}, {1.}, {1.}});
  }

  for(const FixedOrientationContact& fc: fixedOriContacts_)
  {
    boost::shared_ptr<FixedOrientationContactConstr<Type>> fcc(
        new FixedOrientationContactConstr<Type>(&pgdata_, fc.bodyId, fc.target, fc.surfaceFrame));
    problem.addConstraint(fcc, {{0., 0.}, {0., 0.}, {0., 0.}},
        {{1e-2}, {1e-2}, {1e-2}});
  }

  if(!forceContacts_.empty())
  {
    boost::shared_ptr<StaticStabilityConstr<Type>> stab(
        new StaticStabilityConstr<Type>(&pgdata_));
    problem.addConstraint(stab, {{0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}},
        {{1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}});
  }

  roboptim::IpoptSolver solver(problem);

  for(const auto& p: params_)
  {
    solver.parameters()[p.first] = p.second;
  }

  solver_t::result_t res = solver.minimum();
  // Check if the minimization has succeed.
  if(res.which() != solver_t::SOLVER_VALUE &&
     res.which() != solver_t::SOLVER_VALUE_WARNINGS)
  {
    return false;
  }

  ResultVisitor resVisitor;
  boost::apply_visitor(resVisitor, res);
  x_ = resVisitor.x;

  return true;
}


template<typename Type>
std::vector<std::vector<double> > PostureGenerator<Type>::q() const
{
  Eigen::VectorXd eigenQ = x_.head(pgdata_.multibody().nrParams());
  if(pgdata_.multibody().joint(0).type() == rbd::Joint::Free)
  {
    eigenQ.head(4) /= eigenQ.head(4).norm();
  }

  return rbd::vectorToParam(pgdata_.multibody(), eigenQ);
}


} // namespace pg
