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
// boost

// include
// std
#include <fstream>
#include <iostream>
#include <tuple>

// boost
#define BOOST_TEST_MODULE Algo test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// SCD
#include <SCD/S_Object/S_Sphere.h>

#include "EigenAutoDiffScalar.h"
// RBDyn
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

// PG
#include "FK.h"
#include "ID.h"
#include "PGData.h"
#include "PostureGenerator.h"

const Eigen::Vector3d gravity(0., 9.81, 0.);

/// @return An simple ZXZ arm with Y as up axis.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeZXZArm(bool isFixed=true)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  MultiBodyGraph mbg;

  double mass = 1.;
  Matrix3d I = Matrix3d::Identity();
  Vector3d h = Vector3d::Zero();

  RBInertiad rbi(mass, h, I);

  Body b0(rbi, 0, "b0");
  Body b1(rbi, 1, "b1");
  Body b2(rbi, 2, "b2");
  Body b3(rbi, 3, "b3");

  mbg.addBody(b0);
  mbg.addBody(b1);
  mbg.addBody(b2);
  mbg.addBody(b3);

  Joint j0(Joint::RevZ, true, 0, "j0");
  Joint j1(Joint::RevX, true, 1, "j1");
  Joint j2(Joint::RevZ, true, 2, "j2");

  mbg.addJoint(j0);
  mbg.addJoint(j1);
  mbg.addJoint(j2);

  //  Root     j0       j1     j2
  //  ---- b0 ---- b1 ---- b2 ----b3
  //  Fixed    Z       X       Z


  PTransformd to(Vector3d(0., 0.5, 0.));
  PTransformd from(Vector3d(0., 0., 0.));


  mbg.linkBodies(0, PTransformd::Identity(), 1, from, 0);
  mbg.linkBodies(1, to, 2, from, 1);
  mbg.linkBodies(2, to, 3, from, 2);

  MultiBody mb = mbg.makeMultiBody(0, isFixed);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc);
}


/// @return An simple Z*12 arm with Y as up axis.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeZ12Arm(bool isFixed=true)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  MultiBodyGraph mbg;

  double mass = 1.;
  Matrix3d I = Matrix3d::Identity();
  Vector3d h = Vector3d::Zero();

  RBInertiad rbi(mass, h, I);

  for(int i = 0; i < 13; ++i)
  {
    std::stringstream ss;
    ss << "b" << i;
    mbg.addBody({rbi, i, ss.str()});
  }

  for(int i = 0; i < 12; ++i)
  {
    std::stringstream ss;
    ss << "j" << i;
    mbg.addJoint({Joint::RevZ, true, i, ss.str()});
  }

  PTransformd to(Vector3d(0., 0.5, 0.));
  PTransformd from(Vector3d(0., 0., 0.));

  mbg.linkBodies(0, PTransformd::Identity(), 1, from, 0);
  for(int i = 1; i < 12; ++i)
  {
    mbg.linkBodies(i, to, i + 1, from, i);
  }

  MultiBody mb = mbg.makeMultiBody(0, isFixed);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc);
}


BOOST_AUTO_TEST_CASE(FKTest)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  namespace cst = boost::math::constants;
  using scalar_t = pg::eigen_ad::scalar_t;

  MultiBody mb;
  MultiBodyConfig mbcInit;

  std::tie(mb, mbcInit) = makeZXZArm();

  pg::FK<scalar_t> fk(mb);
  std::vector<std::vector<scalar_t>> q = {{},
                                          {scalar_t(0., 3, 0)},
                                          {scalar_t(0., 3, 1)},
                                          {scalar_t(0., 3, 2)}};
  fk.run(mb, q);
}


BOOST_AUTO_TEST_CASE(IDTest)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  namespace cst = boost::math::constants;
  using scalar_t = pg::eigen_ad::scalar_t;

  MultiBody mb;
  MultiBodyConfig mbcInit;

  std::tie(mb, mbcInit) = makeZXZArm();
  double mass = 0.;
  for(const Body& b: mb.bodies())
  {
    mass += b.inertia().mass();
  }

  pg::FK<scalar_t> fk(mb);
  std::vector<std::vector<scalar_t>> q = {{},
                                          {scalar_t(0., 6, 0)},
                                          {scalar_t(0., 6, 1)},
                                          {scalar_t(0., 6, 2)}};
  fk.run(mb, q);

  pg::ID<scalar_t> id(mb, mbcInit.gravity);

  ForceVec<scalar_t> fNull(Vector6<scalar_t>::Zero());
  ForceVec<scalar_t> fbaseNull(Vector3<scalar_t>::Zero(),
                               Vector3<scalar_t>(
                                 scalar_t(0., 6, 3),
                                 scalar_t(mass*9.81, 6, 4),
                                 scalar_t(0., 6, 5)
                                 ));

  id.run(mb, fk.bodyPosW(), fk.parentToSon(), {fbaseNull, fNull, fNull, fNull});
  BOOST_CHECK_SMALL((pg::toValue(id.bodyAcc()[0].vector()) - Vector6d::Zero()).norm(), 1e-6);
}


BOOST_AUTO_TEST_CASE(PGDataTest)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  namespace cst = boost::math::constants;

  MultiBody mb;
  MultiBodyConfig mbcInit;

  std::tie(mb, mbcInit) = makeZXZArm();

  pg::PGData<pg::eigen_ad> pgdata(mb, gravity);
  pgdata.x(Eigen::VectorXd::Zero(3));
}


BOOST_AUTO_TEST_CASE(PGTest)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  namespace cst = boost::math::constants;

  MultiBody mb;
  MultiBodyConfig mbcInit, mbcWork;

  std::tie(mb, mbcInit) = makeZXZArm();
  mbcWork = mbcInit;

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);

    Vector3d target(0., 0.5, 0.5);
    pgPb.fixedPositionContacts({{3, target, sva::PTransformd::Identity()}});

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[3].translation() - target).norm(), 1e-5);
  }

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);

    Matrix3d target(Quaterniond(AngleAxisd(-cst::pi<double>()/2., Vector3d::UnitX())));
    pgPb.fixedOrientationContacts({{3, target, sva::PTransformd::Identity()}});

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[3].rotation() - target).norm(), 1e-3);
  }
}


void toPython(const rbd::MultiBody& mb,
              const rbd::MultiBodyConfig& mbc,
              const std::vector<pg::ForceContact>& fc,
              const std::vector<sva::ForceVecd>& forces,
              const std::string& filename)
{
  std::ofstream out(filename);

  out << "pos = [";
  for(std::size_t i = 0; i < mbc.bodyPosW.size(); ++i)
  {
    out << "[" << mbc.bodyPosW[i].translation()[0] << ", "
               << mbc.bodyPosW[i].translation()[1] << ", "
               << mbc.bodyPosW[i].translation()[2] << "], ";
  }
  out << "]" << std::endl;

  out << "forces = [";
  std::size_t findex = 0;
  for(std::size_t i = 0; i < fc.size(); ++i)
  {
    int bodyIndex = mb.bodyIndexById(fc[i].bodyId);
    for(std::size_t j = 0; j < fc[i].points.size(); ++j)
    {
      Eigen::Vector3d start = (fc[i].points[j]*mbc.bodyPosW[bodyIndex]).translation();
      out << "[";
      out << "(" << start[0] << ", "
                 << start[1] << ", "
                 << start[2] << "), ";
      Eigen::Vector3d end = start + forces[findex].force();
      out << "(" << end[0] << ", "
                 << end[1] << ", "
                 << end[2] << "), ";
      out << "]," << std::endl;
      ++findex;
    }
  }
  out << "]" << std::endl;
}


BOOST_AUTO_TEST_CASE(PGTestZ12)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  namespace cst = boost::math::constants;

  MultiBody mb;
  MultiBodyConfig mbcInit, mbcWork;

  std::tie(mb, mbcInit) = makeZ12Arm();
  // to avoid to start in singularity
  mbcInit.q[3][0] = -0.1;
  mbcWork = mbcInit;

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    Vector3d target(2., 0., 0.);
    Matrix3d oriTarget(sva::RotZ(-cst::pi<double>()));
    int id = 12;
    int index = mb.bodyIndexById(id);
    pgPb.fixedPositionContacts({{id, target, sva::PTransformd::Identity()}});
    pgPb.fixedOrientationContacts({{id, oriTarget, sva::PTransformd::Identity()}});

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].translation() - target).norm(), 1e-5);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].rotation() - oriTarget).norm(), 1e-3);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12.py");
  }

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    Vector3d target(2., 0., 0.);
    Matrix3d oriTarget(sva::RotZ(-cst::pi<double>()));
    int id = 12;
    int index = mb.bodyIndexById(id);
    pgPb.fixedPositionContacts({{id, target, sva::PTransformd::Identity()}});
    pgPb.fixedOrientationContacts({{id, oriTarget, sva::PTransformd::Identity()}});
    Matrix3d frame(RotX(-cst::pi<double>()/2.));
    pgPb.forceContacts({{0, {sva::PTransformd(frame, Vector3d(0.01, 0., 0.)),
                             sva::PTransformd(frame, Vector3d(-0.01, 0., 0.))}, 1.}});

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].translation() - target).norm(), 1e-5);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].rotation() - oriTarget).norm(), 1e-3);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12Stab.py");
  }

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    Vector3d target(1.5, 0., 0.);
    Matrix3d oriTarget(sva::RotZ(-cst::pi<double>()));
    int id = 12;
    int index = mb.bodyIndexById(id);
    pgPb.fixedPositionContacts({{id, target, sva::PTransformd::Identity()}});
    pgPb.fixedOrientationContacts({{id, oriTarget, sva::PTransformd::Identity()}});
    Matrix3d frame(RotX(-cst::pi<double>()/2.));
    Matrix3d frameEnd(RotX(cst::pi<double>()/2.));
    pgPb.forceContacts({{0 , {sva::PTransformd(frame, Vector3d(0.01, 0., 0.)),
                              sva::PTransformd(frame, Vector3d(-0.01, 0., 0.))}, 1.},
                        {id, {sva::PTransformd(frameEnd, Vector3d(0.01, 0., 0.)),
                              sva::PTransformd(frameEnd, Vector3d(-0.01, 0., 0.))}, 1.}});

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].translation() - target).norm(), 1e-5);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].rotation() - oriTarget).norm(), 1e-3);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12Stab2.py");
  }

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    int id = 12;
    int index = mb.bodyIndexById(id);
    Matrix3d frame(RotX(-cst::pi<double>()/2.));
    sva::PTransformd targetSurface(frame, Vector3d(0., 1., 0.));
    sva::PTransformd bodySurface(frame);
    std::vector<Eigen::Vector2d> targetPoints = {{1., 1.}, {-0., 1.}, {-0., -1.}, {1., -1.}};
    std::vector<Eigen::Vector2d> surfPoints = {{0.1, 0.1}, {-0.1, 0.1}, {-0.1, -0.1}, {0.1, -0.1}};
    pgPb.planarContacts({{id, targetSurface, targetPoints, bodySurface, surfPoints}});
    pgPb.bodyPositionTargets({{id, Vector3d(2., 1., 0.), 10.}});

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    sva::PTransformd surfPos = bodySurface*mbcWork.bodyPosW[index];
    double posErr = (surfPos.translation() -
        targetSurface.translation()).dot(targetSurface.rotation().row(2));
    double oriErr = surfPos.rotation().row(2).dot(targetSurface.rotation().row(2)) - 1.;
    BOOST_CHECK_SMALL(posErr, 1e-5);
    BOOST_CHECK_SMALL(oriErr, 1e-5);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12Planar.py");
  }

  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    Vector3d target(1.5, 0., 0.);
    Matrix3d oriTarget(sva::RotZ(-cst::pi<double>()));
    int id = 12;
    pgPb.fixedPositionContacts({{id, target, sva::PTransformd::Identity()}});
    pgPb.fixedOrientationContacts({{id, oriTarget, sva::PTransformd::Identity()}});
    Matrix3d frame(RotX(-cst::pi<double>()/2.));
    Matrix3d frameEnd(RotX(cst::pi<double>()/2.));
    std::vector<pg::ForceContact> fcVec =
        {{0 , {sva::PTransformd(frame, Vector3d(0.01, 0., 0.)),
               sva::PTransformd(frame, Vector3d(-0.01, 0., 0.))}, 1.},
         {id, {sva::PTransformd(frameEnd, Vector3d(0.01, 0., 0.)),
               sva::PTransformd(frameEnd, Vector3d(-0.01, 0., 0.))}, 1.}};
    pgPb.forceContacts(fcVec);

    std::vector<std::vector<double>> ql(mb.nrJoints());
    std::vector<std::vector<double>> qu(mb.nrJoints());
    for(std::size_t i = 0; i < ql.size(); ++i)
    {
      ql[i].resize(mb.joint(int(i)).dof());
      qu[i].resize(mb.joint(int(i)).dof());
      for(std::size_t j = 0; j < ql[i].size(); ++j)
      {
        ql[i][j] = -100.;
        qu[i][j] = 100.;
      }
    }
    pgPb.torqueBounds(ql, qu);

    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    std::vector<sva::ForceVecd> forces = pgPb.forces();
    std::vector<std::vector<double>> torque = pgPb.torque();

    mbcWork.zero(mb);
    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    forwardVelocity(mb, mbcWork);

    // Input force computed by the pg
    int forceIndex = 0;
    for(const pg::ForceContact& f: fcVec)
    {
      int index = mb.bodyIndexById(f.bodyId);
      for(const sva::PTransformd& p: f.points)
      {
        mbcWork.force[index] = mbcWork.force[index] +
            mbcWork.bodyPosW[index].inv().dualMul(p.inv().dualMul(forces[forceIndex]));
        ++forceIndex;
      }
    }

    // Compute the inverse dynamics
    rbd::InverseDynamics invDyn(mb);
    invDyn.inverseDynamics(mb, mbcWork);

    // check if torque are equals
    for(int i = 0; i < mb.nrJoints(); ++i)
    {
      for(int j = 0; j < mb.joint(i).dof(); ++j)
      {
        BOOST_CHECK_SMALL(mbcWork.jointTorque[i][j] -torque[i][j], 1e-5);
      }
    }

    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12Torque.py");
  }


  /*
   *                      Environment collision avoidance
   */
  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    Vector3d target(0., 0., 0.);
    int id = 12;
    int index = mb.bodyIndexById(id);
    pgPb.bodyPositionTargets({{id, target, 0.1}});

    // first we try to go to origin
    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    auto qOrigin = pgPb.q();
    mbcWork.q = qOrigin;
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_SMALL((mbcWork.bodyPosW[index].translation() - target).norm(), 1e-5);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12EnvCol0.py");

    pgPb.param("ipopt.tol", 1e-1);
    pgPb.param("ipopt.dual_inf_tol", 1e-1);
    SCD::S_Sphere hullBody(0.5);
    SCD::S_Sphere hullEnv(0.5);
    hullEnv.setTransformation(pg::toSCD(sva::PTransformd::Identity()));

    pgPb.envCollisions({{id, &hullBody, sva::PTransformd::Identity(), &hullEnv, 0.1}});
    // we check that we couldn't go in collision
    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_GT((mbcWork.bodyPosW[index].translation() - target).norm(), 1. + 0.1);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12EnvCol1.py");


    pgPb.bodyPositionTargets({});
    // same check but we start in constraint violation
    BOOST_REQUIRE(pgPb.run(qOrigin, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    BOOST_CHECK_GT((mbcWork.bodyPosW[index].translation() - target).norm(), 1. + 0.1);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12EnvCol2.py");
  }


  /*
   *                      Self collision avoidance
   */
  {
    pg::PostureGenerator<pg::eigen_ad> pgPb(mb, gravity);
    pgPb.param("ipopt.print_level", 0);
    pgPb.param("ipopt.linear_solver", "ma27");

    Vector3d target(0., 0., 0.);
    int id1 = 12;
    int index1 = mb.bodyIndexById(id1);
    int id2 = 6;
    int index2 = mb.bodyIndexById(id2);
    pgPb.bodyPositionTargets({{id1, target, 0.1}, {id2, target, 0.1}});

    // first we try to go to origin
    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    auto qOrigin = pgPb.q();
    mbcWork.q = qOrigin;
    forwardKinematics(mb, mbcWork);
    double bodyDist = (mbcWork.bodyPosW[index1].translation() -
                       mbcWork.bodyPosW[index2].translation()).norm();
    BOOST_CHECK_SMALL(bodyDist, 1e-5);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12SelfCol0.py");

    pgPb.param("ipopt.tol", 1e-1);
    pgPb.param("ipopt.dual_inf_tol", 1e-1);
    SCD::S_Sphere hullBody1(0.5);
    SCD::S_Sphere hullBody2(0.5);

    pgPb.selfCollisions({{id1, &hullBody1, sva::PTransformd::Identity(),
                          id2, &hullBody2, sva::PTransformd::Identity(),
                          0.1}});

    // we check that we couldn't go in collision
    BOOST_REQUIRE(pgPb.run(mbcInit.q, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    bodyDist = (mbcWork.bodyPosW[index1].translation() -
                mbcWork.bodyPosW[index2].translation()).norm();
    BOOST_CHECK_GT(bodyDist, 1. + 0.1);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12SelfCol1.py");


    pgPb.bodyPositionTargets({});
    // same check but we start in constraint violation
    BOOST_REQUIRE(pgPb.run(qOrigin, {}, mbcInit.q, 0., 0., 0.));

    mbcWork.q = pgPb.q();
    forwardKinematics(mb, mbcWork);
    bodyDist = (mbcWork.bodyPosW[index1].translation() -
                mbcWork.bodyPosW[index2].translation()).norm();
    BOOST_CHECK_GT(bodyDist, 1. + 0.1);
    toPython(mb, mbcWork, pgPb.forceContacts(), pgPb.forces(),"Z12SelfCol2.py");
  }
}
