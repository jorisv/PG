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

namespace pg
{

struct FixedContact
{
  int bodyId;
  Eigen::Vector3d pos;
  /// @todo add friction and contact point in an inner struct member
};


template<typename Type>
class PostureGenerator
{
public:
  PostureGenerator(const rbd::MultiBody& mb);

  void fixedContacts(std::vector<FixedContact> contacts);
  void qBounds(const std::vector<std::vector<double>>& lq,
               const std::vector<std::vector<double>>& uq);

  bool run(const std::vector<std::vector<double>>& q);

  std::vector<std::vector<double>> q() const;

private:
  PGData<Type> pgdata_;

  std::vector<FixedContact> fixedContacts_;
  Eigen::VectorXd ql_, qu_;

  Eigen::VectorXd x_;
};


// inline


template<typename Type>
PostureGenerator<Type>::PostureGenerator(const rbd::MultiBody& mb)
  : pgdata_(mb)
  , fixedContacts_()
  , ql_(mb.nrParams())
  , qu_(mb.nrParams())
{
  ql_.setConstant(-std::numeric_limits<double>::infinity());
  qu_.setConstant(std::numeric_limits<double>::infinity());
}


template<typename Type>
void PostureGenerator<Type>::fixedContacts(std::vector<FixedContact> contacts)
{
  fixedContacts_ = std::move(contacts);
}


template<typename Type>
void PostureGenerator<Type>::qBounds(const std::vector<std::vector<double>>& ql,
    const std::vector<std::vector<double>>& qu)
{
  ql_ = rbd::paramToVector(pgdata_.multibody(), ql);
  qu_ = rbd::paramToVector(pgdata_.multibody(), qu);
}


template<typename Type>
bool PostureGenerator<Type>::run(const std::vector<std::vector<double> >& q)
{
  /*
  typedef roboptim::Solver<roboptim::DifferentiableFunction,
      boost::mpl::vector<roboptim::LinearFunction, roboptim::DifferentiableFunction>> solver_t;
      */
  typedef roboptim::IpoptSolver::solver_t solver_t;

  StdCostFunc<Type> cost(&pgdata_);

  solver_t::problem_t problem(cost);
  problem.startingPoint() = Eigen::VectorXd::Zero(pgdata_.pbSize());
  problem.startingPoint()->head(pgdata_.multibody().nrParams()) =
      rbd::paramToVector(pgdata_.multibody(), q);

  for(int i = 0; i < pgdata_.multibody().nrParams(); ++i)
  {
    problem.argumentBounds()[i] = {ql_[i], qu_[i]};
  }

  for(const FixedContact& fc: fixedContacts_)
  {
    boost::shared_ptr<FixedContactConstr<Type>> fcc(
        new FixedContactConstr<Type>(&pgdata_, fc.bodyId, fc.pos));
    problem.addConstraint(fcc, {{0., 0.}, {0., 0.}, {0., 0.}},
        {{1.}, {1.}, {1.}});
  }

  /*
  roboptim::SolverFactory<solver_t> factory("ipopt", problem);
  solver_t& solver = factory();
  */
  roboptim::IpoptSolver solver(problem);

  solver_t::result_t res = solver.minimum();
  // Check if the minimization has succeed.
  if(res.which () != solver_t::SOLVER_VALUE)
  {
    return false;
  }

  roboptim::Result& result = boost::get<roboptim::Result>(res);
  x_ = result.x;

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
