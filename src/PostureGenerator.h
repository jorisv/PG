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
  typedef roboptim::EigenMatrixDense functionType_t;

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
  PostureGenerator(const rbd::MultiBody& mb, const Eigen::Vector3d& gravity);

  void fixedPositionContacts(std::vector<FixedPositionContact> contacts);
  void fixedOrientationContacts(std::vector<FixedOrientationContact> contacts);
  void planarContacts(std::vector<PlanarContact> contacts);
  void ellipseContacts(std::vector<EllipseContact> contacts);
  void gripperContacts(std::vector<GripperContact> contacts);

  void forceContacts(std::vector<ForceContact> contacts);
  const std::vector<ForceContact>& forceContacts();

  void envCollisions(std::vector<EnvCollision> cols);
  void selfCollisions(std::vector<SelfCollision> cols);

  void bodyPositionTargets(std::vector<BodyPositionTarget> targets);
  void bodyOrientationTargets(std::vector<BodyOrientationTarget> targets);

  void forceContactsMinimization(std::vector<ForceContactMinimization> min);
  void ellipseCostScale(double ellipseScale);

  void qBounds(const std::vector<std::vector<double>>& lq,
               const std::vector<std::vector<double>>& uq);
  void torqueBounds(std::vector<std::vector<double>> lt,
                    std::vector<std::vector<double>> ut);
  void torqueBoundsPoly(std::vector<std::vector<Eigen::VectorXd>> lt,
                        std::vector<std::vector<Eigen::VectorXd>> ut);

  void param(const std::string& name, const std::string& value);
  void param(const std::string& name, double value);
  void param(const std::string& name, int value);

  bool run(const std::vector<std::vector<double> >& initQ,
           const std::vector<sva::ForceVecd>& initForces,
           const std::vector<std::vector<double> >& targetQ,
           double postureScale, double torqueScale, double forceScale);

  std::vector<std::vector<double>> q() const;
  std::vector<sva::ForceVecd> forces() const;
  std::vector<std::vector<double>> torque();
  std::vector<EllipseResult> ellipses() const;

  int nrIters() const;
  std::vector<std::vector<double>> qIter(int i) const;
  std::vector<sva::ForceVecd> forcesIter(int i) const;
  std::vector<std::vector<double>> torqueIter(int i);
  std::vector<EllipseResult> ellipsesIter(int i) const;
  IterateQuantities quantitiesIter(int i) const;

private:
  std::vector<std::vector<double>> q(const Eigen::VectorXd& x) const;
  std::vector<sva::ForceVecd> forces(const Eigen::VectorXd& x) const;
  std::vector<std::vector<double>> torque(const Eigen::VectorXd& x);
  std::vector<EllipseResult> ellipses(const Eigen::VectorXd& x) const;

private:
  PGData pgdata_;
  solver_t::parameters_t params_;

  std::vector<FixedPositionContact> fixedPosContacts_;
  std::vector<FixedOrientationContact> fixedOriContacts_;
  std::vector<PlanarContact> planarContacts_;
  std::vector<EllipseContact> ellipseContacts_;
  double ellipseCostScale_;
  std::vector<GripperContact> gripperContacts_;
  std::vector<ForceContact> forceContacts_;
  std::vector<EnvCollision> envCollisions_;
  std::vector<SelfCollision> selfCollisions_;
  std::vector<BodyPositionTarget> bodyPosTargets_;
  std::vector<BodyOrientationTarget> bodyOriTargets_;
  std::vector<ForceContactMinimization> forceContactsMin_;
  Eigen::VectorXd ql_, qu_;
  Eigen::VectorXd tl_, tu_;
  std::vector<std::vector<Eigen::VectorXd>> tlPoly_, tuPoly_;
  bool isTorque_;

  Eigen::VectorXd x_;
  boost::shared_ptr<iteration_callback_t> iters_;
};


// inline



} // namespace pg
