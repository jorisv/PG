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

// associated header
#include "PostureGenerator.h"

// include
// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// PG
#include "ConfigStruct.h"
#include "PGData.h"
#include "StdCostFunc.h"
#include "FixedContactConstr.h"
#include "StaticStabilityConstr.h"
#include "PositiveForceConstr.h"
#include "FrictionConeConstr.h"
//#include "TorqueConstr.h"
#include "PlanarSurfaceConstr.h"
//#include "EllipseContactConstr.h"
#include "CollisionConstr.h"
#include "IterationCallback.h"

namespace pg
{

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



PostureGenerator::PostureGenerator()
  : pgdata_(nullptr)
  , robotConfig_(nullptr)
  , iters_(new iteration_callback_t)
{}


void PostureGenerator::robotConfig(RobotConfig robotConfig,
  const Eigen::Vector3d& gravity)
{
  robotConfig_.reset(new RobotConfig(robotConfig));
  pgdata_.reset(new PGData(robotConfig_->mb, gravity));
  pgdata_->forces(robotConfig_->forceContacts);
  pgdata_->ellipses(robotConfig_->ellipseContacts);
}


const RobotConfig& PostureGenerator::robotConfig() const
{
  return *robotConfig_;
}


void PostureGenerator::param(const std::string& name, const std::string& value)
{
  params_[name].value = value;
}


void PostureGenerator::param(const std::string& name, double value)
{
  params_[name].value = value;
}


void PostureGenerator::param(const std::string& name, int value)
{
  params_[name].value = value;
}


bool PostureGenerator::run(const std::vector<std::vector<double> >& initQ,
                           const std::vector<sva::ForceVecd>& initForces,
                           const std::vector<std::vector<double> >& targetQ)
{
  StdCostFunc cost(pgdata_.get(), targetQ, robotConfig_->postureScale,
                   robotConfig_->torqueScale,
                   robotConfig_->forceScale,
                   robotConfig_->ellipseCostScale,
                   robotConfig_->bodyPosTargets,
                   robotConfig_->bodyOriTargets,
                   robotConfig_->forceContacts,
                   robotConfig_->forceContactsMin);

  solver_t::problem_t problem(cost);
  problem.startingPoint() = Eigen::VectorXd::Zero(pgdata_->pbSize());
  problem.startingPoint()->head(pgdata_->multibody().nrParams()) =
      rbd::paramToVector(pgdata_->multibody(), initQ);

  // if init force is not well sized we compute it
  if(int(initForces.size()) != pgdata_->nrForcePoints())
  {
    // compute initial force value
    // this is not really smart but work well
    // for contacts with N vector against gravity vector
    double robotMass = 0.;
    for(const rbd::Body& b: pgdata_->multibody().bodies())
    {
      robotMass += b.inertia().mass();
    }

    double initialForce = (pgdata_->gravity().norm()*robotMass)/pgdata_->nrForcePoints();
    int pos = pgdata_->forceParamsBegin();
    for(const PGData::ForceData& fd: pgdata_->forceDatas())
    {
      for(std::size_t i = 0; i < fd.points.size(); ++i)
      {
        sva::PTransformd point = fd.points[i]*pgdata_->mbc().bodyPosW[fd.bodyIndex];
        Eigen::Vector3d force(point.rotation().row(2).transpose()*initialForce);
        (*problem.startingPoint())[pos + 0] = force(0);
        (*problem.startingPoint())[pos + 1] = force(1);
        (*problem.startingPoint())[pos + 2] = force(2);
        pos += 3;
      }
    }
  }
  else
  {
    int pos = pgdata_->forceParamsBegin();
    for(int i = 0; i < pgdata_->nrForcePoints(); ++i)
    {
      (*problem.startingPoint())[pos + 0] = initForces[i].force()[0];
      (*problem.startingPoint())[pos + 1] = initForces[i].force()[1];
      (*problem.startingPoint())[pos + 2] = initForces[i].force()[2];
      pos += 3;
    }
  }

  for(std::size_t i = 0; i < robotConfig_->ql.size(); ++i)
  {
    for(std::size_t j = 0; j < robotConfig_->ql[i].size(); ++j)
    {
      int jointInParam = robotConfig_->mb.jointPosInParam(int(i));
      problem.argumentBounds()[jointInParam + j] = {robotConfig_->ql[i][j],
                                                    robotConfig_->qu[i][j]};
    }
  }

  for(const FixedPositionContact& fc: robotConfig_->fixedPosContacts)
  {
    boost::shared_ptr<FixedPositionContactConstr> fcc(
        new FixedPositionContactConstr(pgdata_.get(), fc.bodyId, fc.target, fc.surfaceFrame));
    problem.addConstraint(fcc, {{0., 0.}, {0., 0.}, {0., 0.}},
        {{1.}, {1.}, {1.}});
  }

  for(const FixedOrientationContact& fc: robotConfig_->fixedOriContacts)
  {
    boost::shared_ptr<FixedOrientationContactConstr> fcc(
        new FixedOrientationContactConstr(pgdata_.get(), fc.bodyId, fc.target, fc.surfaceFrame));
    problem.addConstraint(fcc, {{1., 1.}, {1., 1.}, {1., 1.}},
        {{1e+1}, {1e+1}, {1e+1}});
  }

  for(const PlanarContact& pc: robotConfig_->planarContacts)
  {
    int bodyIndex = pgdata_->mb().bodyIndexById(pc.bodyId);
    // if the root body is a planar contact and the root joint is a planar joint
    // we don't need to add planar position and orientation constraint
    if(bodyIndex != 0 || pgdata_->mb().joint(0).type() != rbd::Joint::Planar)
    {
      boost::shared_ptr<PlanarPositionContactConstr> ppc(
          new PlanarPositionContactConstr(pgdata_.get(), pc.bodyId, pc.targetFrame, pc.surfaceFrame));
      problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

      // N axis must be aligned between target and surface frame.
      boost::shared_ptr<PlanarOrientationContactConstr> poc(
          new PlanarOrientationContactConstr(pgdata_.get(), pc.bodyId,
                                             pc.targetFrame, pc.surfaceFrame,
                                             2));
      problem.addConstraint(poc, {{1., 1.}}, {{1.}});
    }

    boost::shared_ptr<PlanarInclusionConstr> pic(
        new PlanarInclusionConstr(pgdata_.get(), pc.bodyId,
                                        pc.targetFrame, pc.targetPoints,
                                        pc.surfaceFrame, pc.surfacePoints));
    typename PlanarInclusionConstr::intervals_t limInc(
          pic->outputSize(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalInc(pic->outputSize(), 1.);
    problem.addConstraint(pic, limInc, scalInc);
  }

  /*
  for(const EllipseContact& ec: ellipseContacts_)
  {
    int ellipseIndex = 0;
    boost::shared_ptr<PlanarPositionContactConstr> ppc(
        new PlanarPositionContactConstr(pgdata_.get(), ec.bodyId, ec.targetFrame, ec.surfaceFrame));
    problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

    // N axis must be aligned between target and surface frame.
    boost::shared_ptr<PlanarOrientationContactConstr> poc(
        new PlanarOrientationContactConstr(pgdata_.get(), ec.bodyId,
                                                 ec.targetFrame, ec.surfaceFrame,
                                                 2));
    problem.addConstraint(poc, {{1., 1.}}, {{1.}});

    namespace cst = boost::math::constants;
    double inf = std::numeric_limits<double>::infinity();
    int pos = pgdata_->ellipseParamsBegin() + 5*ellipseIndex;
    problem.argumentBounds()[pos + 2] = {-2*cst::pi<double>(),
                                          2*cst::pi<double>()};
    problem.argumentBounds()[pos + 3] = {ec.radiusMin1, inf};
    problem.argumentBounds()[pos + 4] = {ec.radiusMin2, inf};

    ++ellipseIndex;
  }
  */

  /*
  if(!ellipseContacts_.empty())
  {
    boost::shared_ptr<EllipseContactConstr> ecc(
        new EllipseContactConstr(pgdata_.get(), ellipseContacts_));
    typename EllipseContactConstr::intervals_t limInc(
          ecc->outputSize(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalInc(ecc->outputSize(), 1.);
    problem.addConstraint(ecc, limInc, scalInc);
  }
  */

  for(const GripperContact& gc: robotConfig_->gripperContacts)
  {
    boost::shared_ptr<PlanarPositionContactConstr> ppc(
        new PlanarPositionContactConstr(pgdata_.get(), gc.bodyId, gc.targetFrame, gc.surfaceFrame));
    problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

    boost::shared_ptr<FixedOrientationContactConstr> focc(
        new FixedOrientationContactConstr(pgdata_.get(), gc.bodyId,
                                          gc.targetFrame.rotation(), gc.surfaceFrame));
    problem.addConstraint(focc, {{1., 1.}, {1., 1.}, {1., 1.}},
        {{1e+1}, {1e+1}, {1e+1}});

    boost::shared_ptr<PlanarInclusionConstr> pic(
        new PlanarInclusionConstr(pgdata_.get(), gc.bodyId,
                                        gc.targetFrame, gc.targetPoints,
                                        gc.surfaceFrame, gc.surfacePoints));
    typename PlanarInclusionConstr::intervals_t limInc(
          pic->outputSize(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalInc(pic->outputSize(), 1.);
    problem.addConstraint(pic, limInc, scalInc);
  }

  if(!robotConfig_->forceContacts.empty())
  {
    boost::shared_ptr<StaticStabilityConstr> stab(
        new StaticStabilityConstr(pgdata_.get()));
    problem.addConstraint(stab, {{0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}},
        {{1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}});

    boost::shared_ptr<PositiveForceConstr> positiveForce(
        new PositiveForceConstr(pgdata_.get()));
    typename PositiveForceConstr::intervals_t limPositive(
          pgdata_->nrForcePoints(), {0., std::numeric_limits<double>::infinity()});
    typename solver_t::problem_t::scales_t scalPositive(pgdata_->nrForcePoints(), 1.);
    problem.addConstraint(positiveForce, limPositive, scalPositive);

    boost::shared_ptr<FrictionConeConstr> frictionCone(
        new FrictionConeConstr(pgdata_.get()));
    typename FrictionConeConstr::intervals_t limFriction(
          pgdata_->nrForcePoints(), {-std::numeric_limits<double>::infinity(), 0.});
    typename solver_t::problem_t::scales_t scalFriction(pgdata_->nrForcePoints(), 1.);
    problem.addConstraint(frictionCone, limFriction, scalFriction);
  }

  if(!robotConfig_->envCollisions.empty())
  {
    boost::shared_ptr<EnvCollisionConstr> ec(
        new EnvCollisionConstr(pgdata_.get(), robotConfig_->envCollisions));
    typename EnvCollisionConstr::intervals_t limCol(ec->outputSize());
    for(std::size_t i = 0; i < limCol.size(); ++i)
    {
      limCol[i] = {std::pow(robotConfig_->envCollisions[i].minDist, 2),
                   std::numeric_limits<double>::infinity()};
    }
    typename solver_t::problem_t::scales_t scalCol(ec->outputSize(), 1.);
    problem.addConstraint(ec, limCol, scalCol);
  }

  if(!robotConfig_->selfCollisions.empty())
  {
    boost::shared_ptr<SelfCollisionConstr> sc(
        new SelfCollisionConstr(pgdata_.get(), robotConfig_->selfCollisions));
    typename EnvCollisionConstr::intervals_t limCol(sc->outputSize());
    for(std::size_t i = 0; i < limCol.size(); ++i)
    {
      limCol[i] = {std::pow(robotConfig_->selfCollisions[i].minDist, 2),
                   std::numeric_limits<double>::infinity()};
    }
    typename solver_t::problem_t::scales_t scalCol(sc->outputSize(), 1.);
    problem.addConstraint(sc, limCol, scalCol);
  }

  /*
  if(isTorque_)
  {
    // if polynome constraint not set we use the static torque constraint
    if(tlPoly_.size() == 0)
    {
      std::vector<std::vector<double>> tl = robotConfig_->tl;
      std::vector<std::vector<double>> tu = robotConfig_->tu;
      tl[0] = {};
      tu[0] = {};
      Eigen::VectorXd tl_ = rbd::dofToVector(pgdata_->multibody(), tl);
      Eigen::VectorXd tu_ = rbd::dofToVector(pgdata_->multibody(), tu);

      boost::shared_ptr<TorqueConstr> torque(
          new TorqueConstr(pgdata_.get()));
      typename TorqueConstr::intervals_t limTorque(torque->outputSize());
      for(std::size_t i = 0; i < limTorque.size(); ++i)
      {
        limTorque[i] = {tl_[i], tu_[i]};
      }

      typename solver_t::problem_t::scales_t scalTorque(torque->outputSize(), 1.);
      problem.addConstraint(torque, limTorque, scalTorque);
    }
    else
    {
      std::vector<std::vector<Eigen::VectorXd>> tl = robotConfig_->tlPoly;
      std::vector<std::vector<Eigen::VectorXd>> tu = robotConfig_->tuPoly;
      tl[0] = {};
      tu[0] = {};

      boost::shared_ptr<TorquePolyBoundsConstr> torque(
          new TorquePolyBoundsConstr(pgdata_.get(), tl, tu));
      typename TorquePolyBoundsConstr::intervals_t limTorque(torque->outputSize());
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
  */

  roboptim::SolverFactory<solver_t> factory ("ipopt-sparse", problem);
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



std::vector<std::vector<double> > PostureGenerator::q() const
{
  return q(x_);
}



std::vector<sva::ForceVecd> PostureGenerator::forces() const
{
  return forces(x_);
}



std::vector<std::vector<double> > PostureGenerator::torque()
{
  return torque(x_);
}


std::vector<EllipseResult> PostureGenerator::ellipses() const
{
  return ellipses(x_);
}


int PostureGenerator::nrIters() const
{
  return int(iters_->datas.size());
}



std::vector<std::vector<double> > PostureGenerator::qIter(int i) const
{
  return q(iters_->datas.at(i).x);
}



std::vector<sva::ForceVecd> PostureGenerator::forcesIter(int i) const
{
  return forces(iters_->datas.at(i).x);
}



std::vector<std::vector<double> > PostureGenerator::torqueIter(int i)
{
  return torque(iters_->datas.at(i).x);
}


std::vector<EllipseResult> PostureGenerator::ellipsesIter(int i) const
{
  return ellipses(iters_->datas.at(i).x);
}


IterateQuantities PostureGenerator::quantitiesIter(int i) const
{
  const iteration_callback_t::Data& d = iters_->datas.at(i);
  return IterateQuantities{d.obj,  d.constr_viol};
}



std::vector<std::vector<double> > PostureGenerator::q(const Eigen::VectorXd& x) const
{
  Eigen::VectorXd eigenQ(x.head(pgdata_->multibody().nrParams()));
  if(pgdata_->multibody().joint(0).type() == rbd::Joint::Free)
  {
    eigenQ.head(4) /= eigenQ.head(4).norm();
  }

  return rbd::vectorToParam(pgdata_->multibody(), eigenQ);
}



std::vector<sva::ForceVecd> PostureGenerator::forces(const Eigen::VectorXd& x) const
{
  std::vector<sva::ForceVecd> res(pgdata_->nrForcePoints());
  int pos = pgdata_->multibody().nrParams();
  for(int i = 0; i < int(res.size()); ++i)
  {
    res[i] = sva::ForceVecd(Eigen::Vector3d::Zero(), x.segment<3>(pos));
    pos += 3;
  }

  return std::move(res);
}



std::vector<std::vector<double> > PostureGenerator::torque(const Eigen::VectorXd& /* x */)
{
  std::vector<std::vector<double>> res(pgdata_->multibody().nrJoints());
  /*
  pgdata_->x(x);
  const auto& torque  = pgdata_->id().torque();
  std::vector<std::vector<double>> res(pgdata_->multibody().nrJoints());
  for(std::size_t i = 0; i < res.size(); ++i)
  {
    res[i].resize(pgdata_->multibody().joint(int(i)).dof());
    for(std::size_t j = 0; j < res[i].size(); ++j)
    {
      res[i][j] = torque[i](j).value();
    }
  }

  */
  return std::move(res);
}


std::vector<EllipseResult> PostureGenerator::ellipses(const Eigen::VectorXd& x) const
{
  std::vector<EllipseResult> res(pgdata_->ellipseDatas().size());
  int pos = pgdata_->ellipseParamsBegin();
  int ellipseIndex = 0;

  for(const EllipseContact& ec: robotConfig_->ellipseContacts)
  {
    pos = pgdata_->ellipseParamsBegin() + 5*ellipseIndex;
    res[ellipseIndex].bodyIndex = ec.bodyId;
    res[ellipseIndex].x = x[pos];
    res[ellipseIndex].y = x[pos + 1];
    res[ellipseIndex].theta = x[pos + 2];
    res[ellipseIndex].r1 = x[pos + 3];
    res[ellipseIndex].r2 = x[pos + 4];

    ellipseIndex++;
  }

  return std::move(res);
}

} // pg
