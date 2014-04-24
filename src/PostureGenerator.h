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
#include <boost/bind.hpp>

// roboptim
#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-factory.hh>

// PG
#include "ConfigStruct.h"
#include "PGData.h"


namespace pg
{
class MultiBody;
template <class problem_t, class solverState_t>
struct IterationCallback;

class PostureGenerator
{
public:
  typedef roboptim::EigenMatrixSparse functionType_t;

  //TODO: this could be written in a cleaner way
  typedef boost::mpl::push_back<
    boost::mpl::vector< >,
    roboptim::GenericLinearFunction<functionType_t> >::type
    constraints1_t;
  typedef boost::mpl::push_back<
    constraints1_t,
    roboptim::GenericDifferentiableFunction<functionType_t> >::type
    constraints_t;

  //define the solver
  typedef ::roboptim::Solver<
      ::roboptim::GenericDifferentiableFunction<functionType_t>,
      constraints_t
      > solver_t;

  typedef IterationCallback < solver_t::problem_t,
    solver_t::solverState_t > iteration_callback_t;

public:
  PostureGenerator();

  void robotConfig(std::vector<RobotConfig> robotConfigs, const Eigen::Vector3d& gravity);
  const std::vector<RobotConfig>& robotConfig() const;

  void robotLinks(std::vector<RobotLink> robotLinks);
  const std::vector<RobotLink>& robotLinks() const;

  void param(const std::string& name, const std::string& value);
  void param(const std::string& name, double value);
  void param(const std::string& name, int value);

  bool run(const std::vector<RunConfig>& configs);

  // robot 0
  std::vector<std::vector<double>> q() const;
  std::vector<sva::ForceVecd> forces() const;
  std::vector<std::vector<double>> torque();
  std::vector<EllipseResult> ellipses() const;

  std::vector<std::vector<double>> q(int robot) const;
  std::vector<sva::ForceVecd> forces(int robot) const;
  std::vector<std::vector<double>> torque(int robot);
  std::vector<EllipseResult> ellipses(int robot) const;

  int nrIters() const;

  // robot 0
  std::vector<std::vector<double>> qIter(int i) const;
  std::vector<sva::ForceVecd> forcesIter(int i) const;
  std::vector<std::vector<double>> torqueIter(int i);
  std::vector<EllipseResult> ellipsesIter(int i) const;

  std::vector<std::vector<double>> qIter(int robot, int i) const;
  std::vector<sva::ForceVecd> forcesIter(int robot, int i) const;
  std::vector<std::vector<double>> torqueIter(int robot, int i);
  std::vector<EllipseResult> ellipsesIter(int robot, int i) const;

  IterateQuantities quantitiesIter(int i) const;

private:
  std::vector<std::vector<double>> q(int robot, const Eigen::VectorXd& x) const;
  std::vector<sva::ForceVecd> forces(int robot, const Eigen::VectorXd& x) const;
  std::vector<std::vector<double>> torque(int robot, const Eigen::VectorXd& x);
  std::vector<EllipseResult> ellipses(int robot, const Eigen::VectorXd& x) const;

private:
  std::vector<PGData> pgdatas_;
  std::vector<RobotConfig> robotConfigs_;
  std::vector<RobotLink> robotLinks_;
  solver_t::parameters_t params_;

  Eigen::VectorXd x_;
  boost::shared_ptr<iteration_callback_t> iters_;
};


} // namespace pg
