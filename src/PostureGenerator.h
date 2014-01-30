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
#include <boost/math/constants/constants.hpp>
#include <boost/bind.hpp>

// roboptim
#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-factory.hh>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// PG
#include "PGData.h"
#include "StdCostFunc.h"
#include "FixedContactConstr.h"
#include "StaticStabilityConstr.h"
#include "PositiveForceConstr.h"
#include "FrictionConeConstr.h"
#include "TorqueConstr.h"
#include "PlanarSurfaceConstr.h"
#include "EllipseContactConstr.h"
#include "CollisionConstr.h"
#include "ConfigStruct.h"
#include "IterationCallback.h"

namespace pg
{

template<typename Type>
class PostureGenerator
{
public:
  typedef typename Type::scalar_t scalar_t;

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
  PGData<Type> pgdata_;
  solver_t::parameters_t params_;

  std::vector<FixedPositionContact> fixedPosContacts_;
  std::vector<FixedOrientationContact> fixedOriContacts_;
  std::vector<PlanarContact> planarContacts_;
  std::vector<EllipseContact> ellipseContacts_;
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
PostureGenerator<Type>::PostureGenerator(const rbd::MultiBody& mb,
                                         const Eigen::Vector3d& gravity)
  : pgdata_(mb, gravity)
  , ql_(mb.nrParams())
  , qu_(mb.nrParams())
  , tl_(mb.nrDof())
  , tu_(mb.nrDof())
  , isTorque_(false)
  , iters_(new iteration_callback_t)
{
  ql_.setConstant(-std::numeric_limits<double>::infinity());
  qu_.setConstant(std::numeric_limits<double>::infinity());
  tl_.setConstant(-std::numeric_limits<double>::infinity());
  tu_.setConstant(std::numeric_limits<double>::infinity());
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
void PostureGenerator<Type>::planarContacts(std::vector<PlanarContact> contacts)
{
  planarContacts_ = std::move(contacts);
}

template<typename Type>
void PostureGenerator<Type>::ellipseContacts(std::vector<EllipseContact> contacts)
{
  typedef PGData<Type> pgdata_t;
  typedef typename pgdata_t::EllipseData ellipsedata_t;

  ellipseContacts_ = std::move(contacts);
  std::vector<ellipsedata_t> ellipseDatas;
  ellipseDatas.reserve(ellipseContacts_.size());
  for(const EllipseContact& ec: ellipseContacts_)
  {
    ellipseDatas.push_back({pgdata_.multibody().bodyIndexById(ec.bodyId), 0., 0., 0., 0., 0.});
  }
  pgdata_.ellipses(ellipseDatas);
}

template<typename Type>
void PostureGenerator<Type>::gripperContacts(std::vector<GripperContact> contacts)
{
  gripperContacts_ = std::move(contacts);
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
    forceDatas.push_back({pgdata_.multibody().bodyIndexById(fc.bodyId), points, forces, fc.mu});
  }
  pgdata_.forces(forceDatas);
}


template<typename Type>
const std::vector<ForceContact>& PostureGenerator<Type>::forceContacts()
{
  return forceContacts_;
}


template<typename Type>
void PostureGenerator<Type>::envCollisions(std::vector<EnvCollision> cols)
{
  envCollisions_ = std::move(cols);
}


template<typename Type>
void PostureGenerator<Type>::selfCollisions(std::vector<SelfCollision> cols)
{
  selfCollisions_ = std::move(cols);
}



template<typename Type>
void PostureGenerator<Type>::bodyPositionTargets(std::vector<BodyPositionTarget> targets)
{
  bodyPosTargets_ = std::move(targets);
}


template<typename Type>
void PostureGenerator<Type>::bodyOrientationTargets(std::vector<BodyOrientationTarget> targets)
{
  bodyOriTargets_ = std::move(targets);
}


template<typename Type>
void PostureGenerator<Type>::forceContactsMinimization(std::vector<ForceContactMinimization> min)
{
  forceContactsMin_ = std::move(min);
}


template<typename Type>
void PostureGenerator<Type>::qBounds(const std::vector<std::vector<double>>& ql,
    const std::vector<std::vector<double>>& qu)
{
  ql_ = rbd::paramToVector(pgdata_.multibody(), ql);
  qu_ = rbd::paramToVector(pgdata_.multibody(), qu);
}


template<typename Type>
void PostureGenerator<Type>::torqueBounds(std::vector<std::vector<double>> tl,
    std::vector<std::vector<double>> tu)
{
  tl[0] = {};
  tu[0] = {};
  tl_ = rbd::dofToVector(pgdata_.multibody(), tl);
  tu_ = rbd::dofToVector(pgdata_.multibody(), tu);
  isTorque_ = true;
}


template<typename Type>
void PostureGenerator<Type>::torqueBoundsPoly(std::vector<std::vector<Eigen::VectorXd>> tl,
    std::vector<std::vector<Eigen::VectorXd>> tu)
{
  tl[0] = {};
  tu[0] = {};
  tlPoly_ = std::move(tl);
  tuPoly_ = std::move(tu);
  isTorque_ = true;
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
bool PostureGenerator<Type>::run(const std::vector<std::vector<double> >& initQ,
                                 const std::vector<sva::ForceVecd>& initForces,
                                 const std::vector<std::vector<double> >& targetQ,
                                 double postureScale, double torqueScale, double forceScale)
{
  pgdata_.update();

  StdCostFunc<Type> cost(&pgdata_, targetQ, postureScale, torqueScale, forceScale,
                         bodyPosTargets_, bodyOriTargets_,
                         forceContacts_, forceContactsMin_);

  solver_t::problem_t problem(cost);
  problem.startingPoint() = Eigen::VectorXd::Zero(pgdata_.pbSize());
  problem.startingPoint()->head(pgdata_.multibody().nrParams()) =
      rbd::paramToVector(pgdata_.multibody(), initQ);

  // if init force is not well sized we compute it
  if(int(initForces.size()) != pgdata_.nrForcePoints())
  {
    // compute initial force value
    // this is not really smart but work well
    // for contacts with N vector against gravity vector
    double robotMass = 0.;
    for(const rbd::Body& b: pgdata_.multibody().bodies())
    {
      robotMass += b.inertia().mass();
    }

    double initialForce = (pgdata_.gravity().norm()*robotMass)/pgdata_.nrForcePoints();
    int pos = pgdata_.multibody().nrParams();
    for(int i = 0; i < pgdata_.nrForcePoints(); ++i)
    {
      (*problem.startingPoint())[pos + 2] = initialForce;
      pos += 3;
    }
  }
  else
  {
    int pos = pgdata_.multibody().nrParams();
    for(int i = 0; i < pgdata_.nrForcePoints(); ++i)
    {
      (*problem.startingPoint())[pos + 0] = initForces[i].force()[0];
      (*problem.startingPoint())[pos + 1] = initForces[i].force()[1];
      (*problem.startingPoint())[pos + 2] = initForces[i].force()[2];
      pos += 3;
    }
  }

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
    problem.addConstraint(fcc, {{1., 1.}, {1., 1.}, {1., 1.}},
        {{1e+1}, {1e+1}, {1e+1}});
  }

  for(const PlanarContact& pc: planarContacts_)
  {
    boost::shared_ptr<PlanarPositionContactConstr<Type>> ppc(
        new PlanarPositionContactConstr<Type>(&pgdata_, pc.bodyId, pc.targetFrame, pc.surfaceFrame));
    problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

    // N axis must be aligned between target and surface frame.
    boost::shared_ptr<PlanarOrientationContactConstr<Type>> poc(
        new PlanarOrientationContactConstr<Type>(&pgdata_, pc.bodyId,
                                                 pc.targetFrame, pc.surfaceFrame,
                                                 2));
    problem.addConstraint(poc, {{1., 1.}}, {{1.}});

    boost::shared_ptr<PlanarInclusionConstr<Type>> pic(
        new PlanarInclusionConstr<Type>(&pgdata_, pc.bodyId,
                                        pc.targetFrame, pc.targetPoints,
                                        pc.surfaceFrame, pc.surfacePoints));
    typename PlanarInclusionConstr<Type>::intervals_t limInc(
          pic->outputSize(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalInc(pic->outputSize(), 1.);
    problem.addConstraint(pic, limInc, scalInc);
  }

  for(const EllipseContact& ec: ellipseContacts_)
  {
    int ellipseIndex = 0;
    boost::shared_ptr<PlanarPositionContactConstr<Type>> ppc(
        new PlanarPositionContactConstr<Type>(&pgdata_, ec.bodyId, ec.targetFrame, ec.surfaceFrame));
    problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

    // N axis must be aligned between target and surface frame.
    boost::shared_ptr<PlanarOrientationContactConstr<Type>> poc(
        new PlanarOrientationContactConstr<Type>(&pgdata_, ec.bodyId,
                                                 ec.targetFrame, ec.surfaceFrame,
                                                 2));
    problem.addConstraint(poc, {{1., 1.}}, {{1.}});

    namespace cst = boost::math::constants;
    double inf = std::numeric_limits<double>::infinity();
    int pos = pgdata_.ellipseParamsBegin() + 5*ellipseIndex;
    problem.argumentBounds()[pos + 2] = {-2*cst::pi<double>(),
                                          2*cst::pi<double>()};
    problem.argumentBounds()[pos + 3] = {ec.radiusMin, inf};
    problem.argumentBounds()[pos + 4] = {ec.radiusMin, inf};

    ++ellipseIndex;
  }

  if(!ellipseContacts_.empty())
  {
    boost::shared_ptr<EllipseContactConstr<Type>> ecc(
        new EllipseContactConstr<Type>(&pgdata_, ellipseContacts_));
    typename EllipseContactConstr<Type>::intervals_t limInc(
          ecc->outputSize(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalInc(ecc->outputSize(), 1.);
    problem.addConstraint(ecc, limInc, scalInc);
  }

  for(const GripperContact& gc: gripperContacts_)
  {
    boost::shared_ptr<PlanarPositionContactConstr<Type>> ppc(
        new PlanarPositionContactConstr<Type>(&pgdata_, gc.bodyId, gc.targetFrame, gc.surfaceFrame));
    problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

    // N axis must be aligned between target and surface frame.
    boost::shared_ptr<PlanarOrientationContactConstr<Type>> pocN(
        new PlanarOrientationContactConstr<Type>(&pgdata_, gc.bodyId,
                                                 gc.targetFrame, gc.surfaceFrame,
                                                 2));
    problem.addConstraint(pocN, {{1., 1.}}, {{1.}});

    // T axis must be aligned between target and surface frame.
    // (B could be choose also)
    boost::shared_ptr<PlanarOrientationContactConstr<Type>> pocT(
        new PlanarOrientationContactConstr<Type>(&pgdata_, gc.bodyId,
                                                 gc.targetFrame, gc.surfaceFrame,
                                                 0));
    problem.addConstraint(pocT, {{1., 1.}}, {{1.}});

    boost::shared_ptr<PlanarInclusionConstr<Type>> pic(
        new PlanarInclusionConstr<Type>(&pgdata_, gc.bodyId,
                                        gc.targetFrame, gc.targetPoints,
                                        gc.surfaceFrame, gc.surfacePoints));
    typename PlanarInclusionConstr<Type>::intervals_t limInc(
          pic->outputSize(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalInc(pic->outputSize(), 1.);
    problem.addConstraint(pic, limInc, scalInc);
  }

  if(!forceContacts_.empty())
  {
    boost::shared_ptr<StaticStabilityConstr<Type>> stab(
        new StaticStabilityConstr<Type>(&pgdata_));
    problem.addConstraint(stab, {{0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}},
        {{1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}});

    boost::shared_ptr<PositiveForceConstr<Type>> positiveForce(
        new PositiveForceConstr<Type>(&pgdata_));
    typename PositiveForceConstr<Type>::intervals_t limPositive(
          pgdata_.nrForcePoints(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalPositive(pgdata_.nrForcePoints(), 1.);
    problem.addConstraint(positiveForce, limPositive, scalPositive);
    /*
     * constraint seem to converge more quickly than variable bound.
     * maybe scale fault ?
    int pos = pgdata_.multibody().nrParams();
    for(int i = 0; i < pgdata_.nrForcePoints(); ++i)
    {
      double inf = std::numeric_limits<double>::infinity();
      problem.argumentBounds()[pos + 2] = {0., inf};
      pos += 3;
    }
    */

    boost::shared_ptr<FrictionConeConstr<Type>> frictionCone(
        new FrictionConeConstr<Type>(&pgdata_));
    typename FrictionConeConstr<Type>::intervals_t limFriction(
          pgdata_.nrForcePoints(), {-std::numeric_limits<double>::infinity(), 0.});
    typename solver_t::problem_t::scales_t scalFriction(pgdata_.nrForcePoints(), 1.);
    problem.addConstraint(frictionCone, limFriction, scalFriction);
  }

  if(!envCollisions_.empty())
  {
    boost::shared_ptr<EnvCollisionConstr<Type>> ec(
        new EnvCollisionConstr<Type>(&pgdata_, envCollisions_));
    typename EnvCollisionConstr<Type>::intervals_t limCol(ec->outputSize());
    for(std::size_t i = 0; i < limCol.size(); ++i)
    {
      limCol[i] = {envCollisions_[i].minDist, std::numeric_limits<double>::infinity()};
    }
    typename solver_t::problem_t::scales_t scalCol(ec->outputSize(), 1.);
    problem.addConstraint(ec, limCol, scalCol);
  }

  if(!selfCollisions_.empty())
  {
    boost::shared_ptr<SelfCollisionConstr<Type>> sc(
        new SelfCollisionConstr<Type>(&pgdata_, selfCollisions_));
    typename EnvCollisionConstr<Type>::intervals_t limCol(sc->outputSize());
    for(std::size_t i = 0; i < limCol.size(); ++i)
    {
      limCol[i] = {selfCollisions_[i].minDist, std::numeric_limits<double>::infinity()};
    }
    typename solver_t::problem_t::scales_t scalCol(sc->outputSize(), 1.);
    problem.addConstraint(sc, limCol, scalCol);
  }

  if(isTorque_)
  {
    // if polynome constraint not set we use the static torque constraint
    if(tlPoly_.size() == 0)
    {
      boost::shared_ptr<TorqueConstr<Type>> torque(
          new TorqueConstr<Type>(&pgdata_));
      typename TorqueConstr<Type>::intervals_t limTorque(torque->outputSize());
      for(std::size_t i = 0; i < limTorque.size(); ++i)
      {
        limTorque[i] = {tl_[i], tu_[i]};
      }

      typename solver_t::problem_t::scales_t scalTorque(torque->outputSize(), 1.);
      problem.addConstraint(torque, limTorque, scalTorque);
    }
    else
    {
      boost::shared_ptr<TorquePolyBoundsConstr<Type>> torque(
          new TorquePolyBoundsConstr<Type>(&pgdata_, tlPoly_, tuPoly_));
      typename TorquePolyBoundsConstr<Type>::intervals_t limTorque(torque->outputSize());
      for(std::size_t i = 0; i < limTorque.size()/2; ++i)
      {
        // 0 <= torque(q, f) - torqueMin(q)
        limTorque[i] = {0., std::numeric_limits<double>::infinity()};
        // torque(q, f) - torqueMax(q) <= 0
        limTorque[i + limTorque.size()/2] = {-std::numeric_limits<double>::infinity(), 0.};
      }

      typename solver_t::problem_t::scales_t scalTorque(torque->outputSize(), 1.);
      problem.addConstraint(torque, limTorque, scalTorque);
    }
  }

  roboptim::SolverFactory<solver_t> factory ("ipopt", problem);
  solver_t& solver = factory ();

  iters_->datas.clear();
  solver.setIterationCallback(boost::ref(*iters_));

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
  return q(x_);
}


template<typename Type>
std::vector<sva::ForceVecd> PostureGenerator<Type>::forces() const
{
  return forces(x_);
}


template<typename Type>
std::vector<std::vector<double> > PostureGenerator<Type>::torque()
{
  return torque(x_);
}

template<typename Type>
std::vector<EllipseResult> PostureGenerator<Type>::ellipses() const
{
  return ellipses(x_);
}

template<typename Type>
int PostureGenerator<Type>::nrIters() const
{
  return int(iters_->datas.size());
}


template<typename Type>
std::vector<std::vector<double> > PostureGenerator<Type>::qIter(int i) const
{
  return q(iters_->datas.at(i).x);
}


template<typename Type>
std::vector<sva::ForceVecd> PostureGenerator<Type>::forcesIter(int i) const
{
  return forces(iters_->datas.at(i).x);
}


template<typename Type>
std::vector<std::vector<double> > PostureGenerator<Type>::torqueIter(int i)
{
  return torque(iters_->datas.at(i).x);
}

template<typename Type>
std::vector<EllipseResult> PostureGenerator<Type>::ellipsesIter(int i) const
{
  return ellipses(iters_->datas.at(i).x);
}

template<typename Type>
IterateQuantities PostureGenerator<Type>::quantitiesIter(int i) const
{
  const iteration_callback_t::Data& d = iters_->datas.at(i);
  return IterateQuantities{d.obj,  d.constr_viol};
}


template<typename Type>
std::vector<std::vector<double> > PostureGenerator<Type>::q(const Eigen::VectorXd& x) const
{
  Eigen::VectorXd eigenQ(x.head(pgdata_.multibody().nrParams()));
  if(pgdata_.multibody().joint(0).type() == rbd::Joint::Free)
  {
    eigenQ.head(4) /= eigenQ.head(4).norm();
  }

  return rbd::vectorToParam(pgdata_.multibody(), eigenQ);
}


template<typename Type>
std::vector<sva::ForceVecd> PostureGenerator<Type>::forces(const Eigen::VectorXd& x) const
{
  std::vector<sva::ForceVecd> res(pgdata_.nrForcePoints());
  int pos = pgdata_.multibody().nrParams();
  for(int i = 0; i < int(res.size()); ++i)
  {
    res[i] = sva::ForceVecd(Eigen::Vector3d::Zero(), x.segment<3>(pos));
    pos += 3;
  }

  return std::move(res);
}


template<typename Type>
std::vector<std::vector<double> > PostureGenerator<Type>::torque(const Eigen::VectorXd& x)
{
  pgdata_.x(x);
  const auto& torque  = pgdata_.id().torque();
  std::vector<std::vector<double>> res(pgdata_.multibody().nrJoints());
  for(std::size_t i = 0; i < res.size(); ++i)
  {
    res[i].resize(pgdata_.multibody().joint(int(i)).dof());
    for(std::size_t j = 0; j < res[i].size(); ++j)
    {
      res[i][j] = torque[i](j).value();
    }
  }

  return std::move(res);
}

template<typename Type>
std::vector<EllipseResult> PostureGenerator<Type>::ellipses(const Eigen::VectorXd& x) const
{
  std::vector<EllipseResult> res(pgdata_.ellipseDatas().size());
  int pos = pgdata_.ellipseParamsBegin();
  int ellipseIndex = 0;

  for(const EllipseContact& ec: ellipseContacts_)
  {
    pos = pgdata_.ellipseParamsBegin() + 5*ellipseIndex;
    res[ellipseIndex].bodyIndex = pgdata_.multibody().bodyIndexById(ec.bodyId);
    res[ellipseIndex].x = x[pos];
    res[ellipseIndex].y = x[pos + 1];
    res[ellipseIndex].theta = x[pos + 2];
    res[ellipseIndex].r1 = x[pos + 3];
    res[ellipseIndex].r2 = x[pos + 4];

    ellipseIndex++;
  }

  return std::move(res);
}

} // namespace pg
