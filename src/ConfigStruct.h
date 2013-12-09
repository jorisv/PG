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
#include <vector>

// Eigen
#include <Eigen/Core>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>


// forward declaration
namespace SCD
{
class S_Object;
}

namespace pg
{

struct FixedPositionContact
{
  FixedPositionContact() {}
  FixedPositionContact(int bId, const Eigen::Vector3d& t,
                       const sva::PTransformd& sf)
    : bodyId(bId)
    , target(t)
    , surfaceFrame(sf)
  {}

  int bodyId;
  Eigen::Vector3d target; ///< Position target in world coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
};


struct FixedOrientationContact
{
  FixedOrientationContact() {}
  FixedOrientationContact(int bId, const Eigen::Matrix3d& t,
                          const sva::PTransformd& sf)
    : bodyId(bId)
    , target(t)
    , surfaceFrame(sf)
  {}

  int bodyId;
  Eigen::Matrix3d target; ///< Orientation target in world coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
};


struct PlanarContact
{
  PlanarContact() {}
  PlanarContact(int bId,
                const sva::PTransformd& tf, std::vector<Eigen::Vector2d> tp,
                const sva::PTransformd& sf, std::vector<Eigen::Vector2d> sp)
    : bodyId(bId)
    , targetFrame(tf)
    , targetPoints(std::move(tp))
    , surfaceFrame(sf)
    , surfacePoints(std::move(sp))
  {}

  int bodyId;
  sva::PTransformd targetFrame; ///< Target frame in world coordinate.
  std::vector<Eigen::Vector2d> targetPoints; ///< Target surface points in surface coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
  std::vector<Eigen::Vector2d> surfacePoints; ///< Body surface points in surface coordinate.
};


struct EllipseContact
{
  EllipseContact() {}
  EllipseContact(int bId, double rMin,
                 const sva::PTransformd& tf, std::vector<Eigen::Vector2d> tp,
                 const sva::PTransformd& sf, std::vector<Eigen::Vector2d> sp)
    : bodyId(bId)
    , radiusMin(rMin)
    , targetFrame(tf)
    , targetPoints(std::move(tp))
    , surfaceFrame(sf)
    , surfacePoints(std::move(sp))
  {}

  int bodyId;
  double radiusMin;
  sva::PTransformd targetFrame; ///< Target frame in world coordinate.
  std::vector<Eigen::Vector2d> targetPoints; ///< Target surface points in surface coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
  std::vector<Eigen::Vector2d> surfacePoints; ///< Body surface points in surface coordinate.
};


struct GripperContact
{
  GripperContact() {}
  GripperContact(int bId,
                const sva::PTransformd& tf, std::vector<Eigen::Vector2d> tp,
                const sva::PTransformd& sf, std::vector<Eigen::Vector2d> sp)
    : bodyId(bId)
    , targetFrame(tf)
    , targetPoints(std::move(tp))
    , surfaceFrame(sf)
    , surfacePoints(std::move(sp))
  {}

  int bodyId;
  sva::PTransformd targetFrame; ///< Target frame in world coordinate.
  std::vector<Eigen::Vector2d> targetPoints; ///< Target surface points in surface coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
  std::vector<Eigen::Vector2d> surfacePoints; ///< Body surface points in surface coordinate.
};


struct ForceContact
{
  ForceContact() {}
  ForceContact(int bId, std::vector<sva::PTransformd> p, double m)
    : bodyId(bId)
    , points(std::move(p))
    , mu(m)
  {}

  int bodyId;
  std::vector<sva::PTransformd> points;
  double mu;
};


struct EnvCollision
{
  EnvCollision() {}
  EnvCollision(int bId, SCD::S_Object* bHull, const sva::PTransformd& bT,
               SCD::S_Object* eHull,
               double md)
    : bodyId(bId)
    , bodyHull(bHull)
    , bodyT(bT)
    , envHull(eHull)
    , minDist(md)
  {}

  int bodyId;
  SCD::S_Object* bodyHull;
  sva::PTransformd bodyT;
  SCD::S_Object* envHull;
  double minDist;
};


struct SelfCollision
{
  SelfCollision() {}
  SelfCollision(int b1Id, SCD::S_Object* b1Hull, const sva::PTransformd& b1T,
                int b2Id, SCD::S_Object* b2Hull, const sva::PTransformd& b2T,
                double md)
    : body1Id(b1Id)
    , body1Hull(b1Hull)
    , body1T(b1T)
    , body2Id(b2Id)
    , body2Hull(b2Hull)
    , body2T(b2T)
    , minDist(md)
  {}

  int body1Id;
  SCD::S_Object* body1Hull;
  sva::PTransformd body1T;
  int body2Id;
  SCD::S_Object* body2Hull;
  sva::PTransformd body2T;
  double minDist;
};


struct BodyPositionTarget
{
  BodyPositionTarget() {}
  BodyPositionTarget(int bId, const Eigen::Vector3d& t, double s)
    : bodyId(bId)
    , target(t)
    , scale(s)
  {}

  int bodyId;
  Eigen::Vector3d target;
  double scale;
};


struct BodyOrientationTarget
{
  BodyOrientationTarget() {}
  BodyOrientationTarget(int bId, const Eigen::Matrix3d& t, double s)
    : bodyId(bId)
    , target(t)
    , scale(s)
  {}

  int bodyId;
  Eigen::Matrix3d target;
  double scale;
};


struct ForceContactMinimization
{
  ForceContactMinimization() {}
  ForceContactMinimization(int bId, double s)
    : bodyId(bId)
    , scale(s)
  {}

  int bodyId;
  double scale;
};


struct IterateQuantities
{
  double obj, dual_inf, constr_viol, complem, overallError;
};

} // namespace pg
