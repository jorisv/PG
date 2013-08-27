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

namespace pg
{

struct FixedPositionContact
{
  int bodyId;
  Eigen::Vector3d target; ///< Position target in world coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
};


struct FixedOrientationContact
{
  int bodyId;
  Eigen::Matrix3d target; ///< Orientation target in world coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
};


struct PlanarContact
{
  int bodyId;
  sva::PTransformd targetFrame; ///< Target frame in world coordinate.
  std::vector<Eigen::Vector2d> targetPoints; ///< Target surface points in surface coordinate.
  sva::PTransformd surfaceFrame; ///< Body surface frame in body coordinate.
  std::vector<Eigen::Vector2d> surfacePoints; ///< Body surface points in surface coordinate.
};


struct ForceContact
{
  int bodyId;
  std::vector<sva::PTransformd> points;
  double mu;
};


struct BodyPositionTarget
{
  int bodyId;
  Eigen::Vector3d target;
  double scale;
};


struct BodyOrientationTarget
{
  int bodyId;
  Eigen::Matrix3d target;
  double scale;
};

} // namespace pg
