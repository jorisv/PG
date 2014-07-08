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
#include "RobotLinkConstr.h"
#include "CylindricalSurfaceConstr.h"
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
  : pgdatas_()
  , robotConfigs_()
  , iters_(new iteration_callback_t)
{}


void PostureGenerator::robotConfigs(std::vector<RobotConfig> robotConfigs,
  const Eigen::Vector3d& gravity)
{
  robotConfigs_.clear();
  pgdatas_.clear();
  robotConfigs_.reserve(robotConfigs.size());
  pgdatas_.reserve(robotConfigs.size());

  int pbSize = 0;
  std::vector<int> qPos, fPos;
  qPos.reserve(robotConfigs.size() + 1);
  fPos.reserve(robotConfigs.size() + 1);

  qPos.push_back(0);
  fPos.push_back(std::accumulate(robotConfigs.begin(), robotConfigs.end(), 0,
                                 [](int val, const RobotConfig& rc)
                                   {return val + rc.mb.nrParams();}));

  for(std::size_t i = 0; i < robotConfigs.size(); ++i)
  {
    const RobotConfig& rc = robotConfigs[i];
    qPos.push_back(qPos.back() + rc.mb.nrParams());
    int nrForceVar = std::accumulate(rc.forceContacts.begin(), rc.forceContacts.end(), 0,
                                     [](int val, const ForceContact& fc)
                                       {return val + fc.points.size()*3;});
    fPos.push_back(fPos.back() + nrForceVar);
    pbSize += rc.mb.nrParams() + nrForceVar;
  }

  for(std::size_t i = 0; i < robotConfigs.size(); ++i)
  {
    const RobotConfig& rc = robotConfigs[i];
    PGData pgdata(rc.mb, gravity, pbSize, qPos[i], fPos[i]);
    pgdata.forces(rc.forceContacts);
    pgdata.ellipses(rc.ellipseContacts);
    pgdatas_.push_back(pgdata);
  }

  robotConfigs_ = std::move(robotConfigs);
}


const std::vector<RobotConfig>& PostureGenerator::robotConfigs() const
{
  return robotConfigs_;
}


void PostureGenerator::robotLinks(std::vector<RobotLink> robotLinks)
{
  robotLinks_ = std::move(robotLinks);
}


const std::vector<RobotLink>& PostureGenerator::robotLinks() const
{
  return robotLinks_;
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


bool PostureGenerator::run(const std::vector<RunConfig>& configs)
{
  StdCostFunc cost(pgdatas_, robotConfigs_, configs);

  solver_t::problem_t problem(cost);
  problem.startingPoint() = Eigen::VectorXd::Zero(pgdatas_[0].pbSize());

  for(std::size_t robotIndex = 0; robotIndex < robotConfigs_.size(); ++robotIndex)
  {
    const RobotConfig& robotConfig = robotConfigs_[robotIndex];
    const RunConfig& config = configs[robotIndex];
    PGData& pgdata = pgdatas_[robotIndex];

    problem.startingPoint()->segment(pgdata.qParamsBegin(), pgdata.mb().nrParams()) =
      rbd::paramToVector(pgdata.multibody(), config.initQ);
    // run forward kinematics to compute initial force
    pgdata.updateKinematics(config.initQ);

    // if init force is not well sized we compute it
    if(int(config.initForces.size()) != pgdata.nrForcePoints())
    {
      // compute initial force value
      // this is not really smart but work well
      // for contacts with N vector against gravity vector
      double robotMass = pgdata.robotMass();

      double initialForce = (pgdata.gravity().norm()*robotMass)/pgdata.nrForcePoints();
      int pos = pgdata.forceParamsBegin();
      for(const PGData::ForceData& fd: pgdata.forceDatas())
      {
        for(std::size_t i = 0; i < fd.points.size(); ++i)
        {
          sva::PTransformd point = fd.points[i]*pgdata.mbc().bodyPosW[fd.bodyIndex];
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
      int pos = pgdata.forceParamsBegin();
      for(int i = 0; i < pgdata.nrForcePoints(); ++i)
      {
        (*problem.startingPoint())[pos + 0] = config.initForces[i].force()[0];
        (*problem.startingPoint())[pos + 1] = config.initForces[i].force()[1];
        (*problem.startingPoint())[pos + 2] = config.initForces[i].force()[2];
        pos += 3;
      }
    }


    for(std::size_t i = 0; i < robotConfig.ql.size(); ++i)
    {
      for(std::size_t j = 0; j < robotConfig.ql[i].size(); ++j)
      {
        int jointInParam = robotConfig.mb.jointPosInParam(int(i));
        problem.argumentBounds()[pgdata.qParamsBegin() + jointInParam + j] =
          {robotConfig.ql[i][j], robotConfig.qu[i][j]};
      }
    }

    for(const FixedPositionContact& fc: robotConfig.fixedPosContacts)
    {
      int bodyIndex = pgdata.mb().bodyIndexById(fc.bodyId);
      // if the root body is a fixed contact and the root joint is fixed
      // this constraint is useless
      if(bodyIndex != 0 || pgdata.mb().joint(0).type() != rbd::Joint::Fixed)
      {
        boost::shared_ptr<FixedPositionContactConstr> fcc(
            new FixedPositionContactConstr(&pgdata, fc.bodyId, fc.target, fc.surfaceFrame));
        problem.addConstraint(fcc, {{0., 0.}, {0., 0.}, {0., 0.}},
            {{1.}, {1.}, {1.}});
      }
    }

    for(const FixedOrientationContact& fc: robotConfig.fixedOriContacts)
    {
      int bodyIndex = pgdata.mb().bodyIndexById(fc.bodyId);
      // if the root body is a fixed contact and the root joint is fixed
      // this constraint is useless
      if(bodyIndex != 0 || pgdata.mb().joint(0).type() != rbd::Joint::Fixed)
      {
        boost::shared_ptr<FixedOrientationContactConstr> fcc(
            new FixedOrientationContactConstr(&pgdata, fc.bodyId, fc.target, fc.surfaceFrame));
        problem.addConstraint(fcc, {{1., 1.}, {1., 1.}, {1., 1.}},
            {{1e+1}, {1e+1}, {1e+1}});
      }
    }

    for(const PlanarContact& pc: robotConfig.planarContacts)
    {
      int bodyIndex = pgdata.mb().bodyIndexById(pc.bodyId);
      // if the root body is a planar contact and the root joint is a planar joint
      // we don't need to add planar position and orientation constraint
      if(bodyIndex != 0 || pgdata.mb().joint(0).type() != rbd::Joint::Planar)
      {
        boost::shared_ptr<PlanarPositionContactConstr> ppc(
            new PlanarPositionContactConstr(&pgdata, pc.bodyId, pc.targetFrame, pc.surfaceFrame));
        problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

        // N axis must be aligned between target and surface frame.
        boost::shared_ptr<PlanarOrientationContactConstr> poc(
            new PlanarOrientationContactConstr(&pgdata, pc.bodyId,
                                               pc.targetFrame, pc.surfaceFrame,
                                               2));
        problem.addConstraint(poc, {{0., std::numeric_limits<double>::infinity()},
                                    {0., 0.}, {0., 0.}},
                                   {{1.}, {1.}, {1.}});
      }

      // add planar inclusion only if there is points in both surfaces
      if(pc.targetPoints.size() != 0 && pc.surfacePoints.size() != 0)
      {
         boost::shared_ptr<PlanarInclusionConstr> pic(
            new PlanarInclusionConstr(&pgdata, pc.bodyId,
                                            pc.targetFrame, pc.targetPoints,
                                            pc.surfaceFrame, pc.surfacePoints));
        typename PlanarInclusionConstr::intervals_t limInc(
              pic->outputSize(), {0., std::numeric_limits<double>::infinity()});
        typename solver_t::problem_t::scales_t scalInc(pic->outputSize(), 1.);
        problem.addConstraint(pic, limInc, scalInc);
      }
    }

    /*
    for(const EllipseContact& ec: ellipseContacts_)
    {
      int ellipseIndex = 0;
      boost::shared_ptr<PlanarPositionContactConstr> ppc(
          new PlanarPositionContactConstr(&pgdata, ec.bodyId, ec.targetFrame, ec.surfaceFrame));
      problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

      // N axis must be aligned between target and surface frame.
      boost::shared_ptr<PlanarOrientationContactConstr> poc(
          new PlanarOrientationContactConstr(&pgdata, ec.bodyId,
                                                   ec.targetFrame, ec.surfaceFrame,
                                                   2));
      problem.addConstraint(poc, {{1., 1.}}, {{1.}});

      namespace cst = boost::math::constants;
      double inf = std::numeric_limits<double>::infinity();
      int pos = pgdata.ellipseParamsBegin() + 5*ellipseIndex;
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
          new EllipseContactConstr(&pgdata, ellipseContacts_));
      typename EllipseContactConstr::intervals_t limInc(
            ecc->outputSize(), {0., std::numeric_limits<double>::infinity()});
      typename solver_t::problem_t::scales_t scalInc(ecc->outputSize(), 1.);
      problem.addConstraint(ecc, limInc, scalInc);
    }
    */

    for(const GripperContact& gc: robotConfig.gripperContacts)
    {
      boost::shared_ptr<PlanarPositionContactConstr> ppc(
          new PlanarPositionContactConstr(&pgdata, gc.bodyId, gc.targetFrame, gc.surfaceFrame));
      problem.addConstraint(ppc, {{0., 0.}}, {{1.}});

      boost::shared_ptr<FixedOrientationContactConstr> focc(
          new FixedOrientationContactConstr(&pgdata, gc.bodyId,
                                            gc.targetFrame.rotation(), gc.surfaceFrame));
      problem.addConstraint(focc, {{1., 1.}, {1., 1.}, {1., 1.}},
          {{1e+1}, {1e+1}, {1e+1}});

      boost::shared_ptr<PlanarInclusionConstr> pic(
          new PlanarInclusionConstr(&pgdata, gc.bodyId,
                                          gc.targetFrame, gc.targetPoints,
                                          gc.surfaceFrame, gc.surfacePoints));
      typename PlanarInclusionConstr::intervals_t limInc(
            pic->outputSize(), {0., std::numeric_limits<double>::infinity()});
      typename solver_t::problem_t::scales_t scalInc(pic->outputSize(), 1.);
      problem.addConstraint(pic, limInc, scalInc);
    }

    for(const CylindricalContact& pc: robotConfig.cylindricalContacts)
    {
      int bodyIndex = pgdata.mb().bodyIndexById(pc.bodyId);
      // if the root body is a cylindrical contact and the root joint is a
      // cylindrical joint we don't need to add cylindrical position and
      // orientation constraint
      if(bodyIndex != 0 || pgdata.mb().joint(0).type() != rbd::Joint::Cylindrical)
      {
        sva::PTransformd radiusX(Eigen::Vector3d(0., 0., pc.targetRadius));
        boost::shared_ptr<CylindricalPositionConstr> fgpc(
            new CylindricalPositionConstr(&pgdata, pc.bodyId, pc.targetFrame,
                                          radiusX*pc.surfaceFrame));
        problem.addConstraint(fgpc, {{-pc.targetWidth/2., pc.targetWidth/2.},
                                     {0., 0.}, {0., 0.}},
                                     {{1.}, {1.}, {1.}});

        // T axis must be aligned between target and surface frame.
        boost::shared_ptr<PlanarOrientationContactConstr> poc(
            new PlanarOrientationContactConstr(&pgdata, pc.bodyId,
                                               pc.targetFrame, pc.surfaceFrame,
                                               0));
        problem.addConstraint(poc, {{0., std::numeric_limits<double>::infinity()},
                                    {0., 0.}, {0., 0.}},
                                    {{1.}, {1.}, {1.}});
      }
    }

    if(!robotConfig.forceContacts.empty())
    {
      boost::shared_ptr<StaticStabilityConstr> stab(
          new StaticStabilityConstr(&pgdata));
      problem.addConstraint(stab, {{0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}},
          {{1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}, {1e-2}});

      boost::shared_ptr<PositiveForceConstr> positiveForce(
          new PositiveForceConstr(&pgdata));
      typename PositiveForceConstr::intervals_t limPositive(
            pgdata.nrForcePoints(), {0., std::numeric_limits<double>::infinity()});
      typename solver_t::problem_t::scales_t scalPositive(pgdata.nrForcePoints(), 1.);
      problem.addConstraint(positiveForce, limPositive, scalPositive);

      boost::shared_ptr<FrictionConeConstr> frictionCone(
          new FrictionConeConstr(&pgdata));
      typename FrictionConeConstr::intervals_t limFriction(
            pgdata.nrForcePoints(), {-std::numeric_limits<double>::infinity(), 0.});
      typename solver_t::problem_t::scales_t scalFriction(pgdata.nrForcePoints(), 1.);
      problem.addConstraint(frictionCone, limFriction, scalFriction);
    }

    if(!robotConfig.envCollisions.empty())
    {
      boost::shared_ptr<EnvCollisionConstr> ec(
          new EnvCollisionConstr(&pgdata, robotConfig.envCollisions));
      typename EnvCollisionConstr::intervals_t limCol(ec->outputSize());
      for(std::size_t i = 0; i < limCol.size(); ++i)
      {
        limCol[i] = {std::pow(robotConfig.envCollisions[i].minDist, 2),
                     std::numeric_limits<double>::infinity()};
      }
      typename solver_t::problem_t::scales_t scalCol(ec->outputSize(), 1.);
      problem.addConstraint(ec, limCol, scalCol);
    }

    if(!robotConfig.selfCollisions.empty())
    {
      boost::shared_ptr<SelfCollisionConstr> sc(
          new SelfCollisionConstr(&pgdata, robotConfig.selfCollisions));
      typename EnvCollisionConstr::intervals_t limCol(sc->outputSize());
      for(std::size_t i = 0; i < limCol.size(); ++i)
      {
        limCol[i] = {std::pow(robotConfig.selfCollisions[i].minDist, 2),
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
        std::vector<std::vector<double>> tl = robotConfig.tl;
        std::vector<std::vector<double>> tu = robotConfig.tu;
        tl[0] = {};
        tu[0] = {};
        Eigen::VectorXd tl_ = rbd::dofToVector(pgdata.multibody(), tl);
        Eigen::VectorXd tu_ = rbd::dofToVector(pgdata.multibody(), tu);

        boost::shared_ptr<TorqueConstr> torque(
            new TorqueConstr(&pgdata));
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
        std::vector<std::vector<Eigen::VectorXd>> tl = robotConfig.tlPoly;
        std::vector<std::vector<Eigen::VectorXd>> tu = robotConfig.tuPoly;
        tl[0] = {};
        tu[0] = {};

        boost::shared_ptr<TorquePolyBoundsConstr> torque(
            new TorquePolyBoundsConstr(&pgdata, tl, tu));
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

    // force PGData to make an update
    // this avoid that a first call with a x identical to PGData::{xq_,xf_}
    // don't update PGData
    pgdata.update();
  }

  for(const RobotLink& rl: robotLinks_)
  {
    boost::shared_ptr<RobotLinkConstr> rlc(
      new RobotLinkConstr(&pgdatas_[rl.robot1Index],
                          &pgdatas_[rl.robot2Index], rl.linkedBodies));
    typename RobotLinkConstr::intervals_t interval(rlc->outputSize(), {0., 0.});
    typename solver_t::problem_t::scales_t scale(rlc->outputSize(), 1.);
    for(std::size_t i = 0; i < rl.linkedBodies.size(); ++i)
    {
      interval[i*6 + 0] = {1., 1.};
      interval[i*6 + 1] = {1., 1.};
      interval[i*6 + 2] = {1., 1.};
      scale[i*6 + 0] = 1e+1;
      scale[i*6 + 1] = 1e+1;
      scale[i*6 + 2] = 1e+1;
    }

    problem.addConstraint(rlc, interval, scale);
  }

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
  return q(0, x_);
}


std::vector<sva::ForceVecd> PostureGenerator::forces() const
{
  return forces(0, x_);
}


std::vector<std::vector<double> > PostureGenerator::torque()
{
  return torque(0, x_);
}


std::vector<EllipseResult> PostureGenerator::ellipses() const
{
  return ellipses(0, x_);
}


std::vector<std::vector<double> > PostureGenerator::q(int robot) const
{
  return q(robot, x_);
}


std::vector<sva::ForceVecd> PostureGenerator::forces(int robot) const
{
  return forces(robot, x_);
}


std::vector<std::vector<double> > PostureGenerator::torque(int robot)
{
  return torque(robot, x_);
}


std::vector<EllipseResult> PostureGenerator::ellipses(int robot) const
{
  return ellipses(robot, x_);
}


int PostureGenerator::nrIters() const
{
  return int(iters_->datas.size());
}


std::vector<std::vector<double> > PostureGenerator::qIter(int i) const
{
  return q(0, iters_->datas.at(i).x);
}



std::vector<sva::ForceVecd> PostureGenerator::forcesIter(int i) const
{
  return forces(0, iters_->datas.at(i).x);
}



std::vector<std::vector<double> > PostureGenerator::torqueIter(int i)
{
  return torque(0, iters_->datas.at(i).x);
}


std::vector<EllipseResult> PostureGenerator::ellipsesIter(int i) const
{
  return ellipses(0, iters_->datas.at(i).x);
}


std::vector<std::vector<double> > PostureGenerator::qIter(int robot, int i) const
{
  return q(robot, iters_->datas.at(i).x);
}


std::vector<sva::ForceVecd> PostureGenerator::forcesIter(int robot, int i) const
{
  return forces(robot, iters_->datas.at(i).x);
}


std::vector<std::vector<double> > PostureGenerator::torqueIter(int robot, int i)
{
  return torque(robot, iters_->datas.at(i).x);
}


std::vector<EllipseResult> PostureGenerator::ellipsesIter(int robot, int i) const
{
  return ellipses(robot, iters_->datas.at(i).x);
}


IterateQuantities PostureGenerator::quantitiesIter(int i) const
{
  const iteration_callback_t::Data& d = iters_->datas.at(i);
  return IterateQuantities{d.obj,  d.constr_viol};
}


std::vector<std::vector<double> >
PostureGenerator::q(int robot, const Eigen::VectorXd& x) const
{
  const PGData& pgdata = pgdatas_[robot];
  Eigen::VectorXd eigenQ(x.segment(pgdata.qParamsBegin(), pgdata.mb().nrParams()));
  if(pgdata.multibody().joint(0).type() == rbd::Joint::Free)
  {
    eigenQ.head(4) /= eigenQ.head(4).norm();
  }

  return rbd::vectorToParam(pgdata.mb(), eigenQ);
}


std::vector<sva::ForceVecd>
PostureGenerator::forces(int robot, const Eigen::VectorXd& x) const
{
  const PGData& pgdata = pgdatas_[robot];
  std::vector<sva::ForceVecd> res(pgdata.nrForcePoints());
  int pos = pgdata.forceParamsBegin();
  for(int i = 0; i < int(res.size()); ++i)
  {
    res[i] = sva::ForceVecd(Eigen::Vector3d::Zero(), x.segment<3>(pos));
    pos += 3;
  }

  return std::move(res);
}


std::vector<std::vector<double> >
PostureGenerator::torque(int robot, const Eigen::VectorXd& /* x */)
{
  const PGData& pgdata = pgdatas_[robot];
  std::vector<std::vector<double>> res(pgdata.multibody().nrJoints());
  /*
  pgdata.x(x);
  const auto& torque  = pgdata.id().torque();
  std::vector<std::vector<double>> res(pgdata.multibody().nrJoints());
  for(std::size_t i = 0; i < res.size(); ++i)
  {
    res[i].resize(pgdata.multibody().joint(int(i)).dof());
    for(std::size_t j = 0; j < res[i].size(); ++j)
    {
      res[i][j] = torque[i](j).value();
    }
  }

  */
  return std::move(res);
}


std::vector<EllipseResult>
PostureGenerator::ellipses(int robot, const Eigen::VectorXd& x) const
{
  const PGData& pgdata = pgdatas_[robot];
  const RobotConfig& robotConfig = robotConfigs_[robot];
  std::vector<EllipseResult> res(pgdata.ellipseDatas().size());
  int pos = pgdata.ellipseParamsBegin();
  int ellipseIndex = 0;

  for(const EllipseContact& ec: robotConfig.ellipseContacts)
  {
    pos = pgdata.ellipseParamsBegin() + 5*ellipseIndex;
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
