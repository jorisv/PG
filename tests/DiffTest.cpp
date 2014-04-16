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

// include
// std
#include <tuple>

// boost
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Diff test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// roboptim
#include <roboptim/core/finite-difference-gradient.hh>

// SCD
#include <SCD/S_Object/S_Sphere.h>

// PG
#include "ConfigStruct.h"
#include "PGData.h"
#include "FixedContactConstr.h"
#include "PlanarSurfaceConstr.h"
#include "StaticStabilityConstr.h"
#include "PositiveForceConstr.h"
#include "FrictionConeConstr.h"
#include "CollisionConstr.h"
#include "StdCostFunc.h"

// Arm
#include "XYZ12Arm.h"

const Eigen::Vector3d gravity(0., 9.81, 0.);

template <typename T>
double checkGradient(
  const roboptim::GenericDifferentiableFunction<T>& function,
  const typename roboptim::GenericDifferentiableFunction<T>::vector_t& x,
  double epsilon=1e-8)
{
  typename roboptim::GenericFunctionTraits<T>::jacobian_t jac(
    std::get<0>(function.jacobianSize()),
    std::get<1>(function.jacobianSize()));
  typename roboptim::GenericFunctionTraits<T>::jacobian_t jacF(
    std::get<0>(function.jacobianSize()),
    std::get<1>(function.jacobianSize()));

  roboptim::GenericFiniteDifferenceGradient<T> fdfunction(function, epsilon);
  function.jacobian(jac, x);
  fdfunction.jacobian(jacF, x);

  return (jac - jacF).norm();
}


template <typename T>
double checkForceGradient(
 const roboptim::GenericDifferentiableFunction<T>& function,
 const typename roboptim::GenericDifferentiableFunction<T>::vector_t& x,
    const pg::PGData& pgdata)
{
  auto rows = std::get<0>(function.jacobianSize());
  auto cols = std::get<1>(function.jacobianSize());
  Eigen::MatrixXd jac(rows, cols);
  Eigen::MatrixXd jacF(rows, cols);

  roboptim::GenericFiniteDifferenceGradient<T> fdfunction(function);
  function.jacobian(jac, x);
  fdfunction.jacobian(jacF, x);

  int deb = pgdata.forceParamsBegin();
  return (jac.block(0, deb, rows, cols - deb) -
          jacF.block(0, deb, rows, cols - deb)).norm();
}


BOOST_AUTO_TEST_CASE(FixedContactPosTest)
{
  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  Eigen::Vector3d target(2., 0., 0.);
  sva::PTransformd surface(sva::PTransformd::Identity());

  pg::FixedPositionContactConstr fpc(&pgdata, 12, target, surface);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(mb.nrDof()));
    BOOST_CHECK_SMALL(checkGradient(fpc, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(FixedContactOriTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  Matrix3d oriTarget(sva::RotZ(cst::pi<double>()));
  sva::PTransformd surface(sva::RotZ(-cst::pi<double>()/2.), Vector3d::Random());

  pg::FixedOrientationContactConstr foc(&pgdata, 3, oriTarget, surface);

  for(int i = 0; i < 100; ++i)
  {
    VectorXd x(VectorXd::Random(mb.nrDof()));
    BOOST_CHECK_SMALL(checkGradient(foc, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(PlanarPositionContactTest)
{
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  sva::PTransformd target(Eigen::Vector3d(0., 1., 0.));
  sva::PTransformd surface(sva::PTransformd::Identity());

  pg::PlanarPositionContactConstr ppp(&pgdata, 12, target, surface);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(mb.nrDof()));
    BOOST_CHECK_SMALL(checkGradient(ppp, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(PlanarOrientationContactTest)
{
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  Eigen::Matrix3d oriTarget(sva::RotZ(-cst::pi<double>()));
  sva::PTransformd surface(sva::PTransformd::Identity());

  pg::PlanarOrientationContactConstr pop(&pgdata, 12, oriTarget, surface, 1);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(mb.nrDof()));
    BOOST_CHECK_SMALL(checkGradient(pop, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(PlanarInclusionTest)
{
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  sva::PTransformd targetSurface(sva::RotZ(-cst::pi<double>()), Eigen::Vector3d(0., 1., 0.));
  sva::PTransformd bodySurface(sva::RotZ(-cst::pi<double>()/2.), Eigen::Vector3d(0., 1., 0.));
  std::vector<Eigen::Vector2d> targetPoints = {{1., 1.}, {-0., 1.}, {-0., -1.}, {1., -1.}};
  std::vector<Eigen::Vector2d> surfPoints = {{0.1, 0.1}, {-0.1, 0.1}, {-0.1, -0.1}, {0.1, -0.1}};

  pg::PlanarInclusionConstr pi(&pgdata, 12, targetSurface, targetPoints,
                               bodySurface, surfPoints);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(mb.nrDof()));
    BOOST_CHECK_SMALL(checkGradient(pi, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(StaticStabilityTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  sva::PTransformd bodySurface(sva::RotZ(-cst::pi<double>()/2.), Eigen::Vector3d(0., 1., 0.));
  std::vector<Vector2d> surfPoints = {{0.1, 0.1}, {-0.1, 0.1}, {-0.1, -0.1}, {0.1, -0.1}};
  std::vector<sva::PTransformd> points(surfPoints.size());
  for(std::size_t i = 0; i < points.size(); ++i)
  {
    points[i] = sva::PTransformd(Vector3d(surfPoints[i][0], surfPoints[i][1], 0.))*bodySurface;
  }
  pgdata.forces({pg::ForceContact{12, points, 0.7}, pg::ForceContact{0, points, 0.7}});

  pg::StaticStabilityConstr ss(&pgdata);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(pgdata.pbSize()));
    BOOST_CHECK_SMALL(checkGradient(ss, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(PositiveForceTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  sva::PTransformd bodySurface(sva::RotZ(-cst::pi<double>()/2.), Eigen::Vector3d(0., 1., 0.));
  std::vector<Vector2d> surfPoints = {{0.1, 0.1}, {-0.1, 0.1}, {-0.1, -0.1}, {0.1, -0.1}};
  std::vector<sva::PTransformd> points(surfPoints.size());
  for(std::size_t i = 0; i < points.size(); ++i)
  {
    points[i] = sva::PTransformd(Vector3d(surfPoints[i][0], surfPoints[i][1], 0.))*bodySurface;
  }
  pgdata.forces({pg::ForceContact{12, points, 0.7}, pg::ForceContact{0, points, 0.7}});

  pg::PositiveForceConstr ss(&pgdata);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(pgdata.pbSize()));
    BOOST_CHECK_SMALL(checkGradient(ss, x), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(FrictionConeTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);

  sva::PTransformd bodySurface(sva::RotZ(-cst::pi<double>()/2.), Eigen::Vector3d(0., 1., 0.));
  std::vector<Vector2d> surfPoints = {{0.1, 0.1}, {-0.1, 0.1}, {-0.1, -0.1}, {0.1, -0.1}};
  std::vector<sva::PTransformd> points(surfPoints.size());
  for(std::size_t i = 0; i < points.size(); ++i)
  {
    points[i] = sva::PTransformd(Vector3d(surfPoints[i][0], surfPoints[i][1], 0.))*bodySurface;
  }
  pgdata.forces({pg::ForceContact{12, points, 0.7}, pg::ForceContact{0, points, 0.7}});

  pg::FrictionConeConstr fc(&pgdata);

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(pgdata.pbSize()));
    BOOST_CHECK_SMALL(checkGradient(fc, x), 1e-4);
  }
}


BOOST_AUTO_TEST_CASE(EnvCollisionTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);
  SCD::S_Sphere hullBody(0.5);
  SCD::S_Sphere hullEnv(0.5);
  hullEnv.setTransformation(pg::toSCD(sva::PTransformd::Identity()));
  pg::EnvCollision ec(12, &hullBody, sva::PTransformd::Identity(),
                      &hullEnv, 0.1);

  pg::EnvCollisionConstr ecc(&pgdata, {ec});

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(pgdata.pbSize())*3.14);
    BOOST_CHECK_SMALL(checkGradient(ecc, x, 1e-4), 1e-1);
  }
}


BOOST_AUTO_TEST_CASE(SelfCollisionTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);
  SCD::S_Sphere hullBody1(0.2);
  SCD::S_Sphere hullBody2(0.2);
  pg::SelfCollision sc(12, &hullBody1, sva::PTransformd::Identity(),
                       6, &hullBody2, sva::PTransformd::Identity(),
                       0.1);

  pg::SelfCollisionConstr scc(&pgdata, {sc});

  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(pgdata.pbSize())*3.14);
    BOOST_CHECK_SMALL(checkGradient(scc, x, 1e-4), 1e-1);
  }
}


BOOST_AUTO_TEST_CASE(StdCostFunctionTest)
{
  using namespace Eigen;
  namespace cst = boost::math::constants;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::tie(mb, mbc) = makeXYZ12Arm();

  pg::PGData pgdata(mb, gravity);
  Vector3d target(1.5, 0., 0.);
  Matrix3d oriTarget(sva::RotZ(-cst::pi<double>()));
  int id = 12;
  Matrix3d frame(sva::RotX(-cst::pi<double>()/2.));
  Matrix3d frameEnd(sva::RotX(cst::pi<double>()/2.));
  std::vector<pg::ForceContact> forceContacts =
    {{0 , {sva::PTransformd(frame, Vector3d(0.01, 0., 0.)),
           sva::PTransformd(frame, Vector3d(-0.01, 0., 0.))}, 1.},
     {id, {sva::PTransformd(frameEnd, Vector3d(0.01, 0., 0.)),
      sva::PTransformd(frameEnd, Vector3d(-0.01, 0., 0.))}, 1.}};
  pgdata.forces(forceContacts);

  std::vector<std::vector<double>> tq(mbc.q);
  for(int i = 1; i < mb.nrJoints(); ++i)
  {
    tq[i] = {0.33*i};
  }

  std::vector<pg::BodyPositionTarget> bodiesPos = {{12, target, 3.45}};
  std::vector<pg::BodyOrientationTarget> bodiesOri = {{12, oriTarget, 5.06}};
  std::vector<pg::ForceContactMinimization> forceMin = {{12, 1.87}};

  pg::StdCostFunc cost(&pgdata, tq, 1.44, 0., 2.33, 0., bodiesPos, bodiesOri,
                       forceContacts, forceMin);


  for(int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd x(Eigen::VectorXd::Random(pgdata.pbSize())*3.14);
    BOOST_CHECK_SMALL(checkGradient(cost, x, 1e-6), 1e-3);
  }
}
