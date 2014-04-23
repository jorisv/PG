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
#include <cassert>
#include <vector>

// boost
#include <boost/math/constants/constants.hpp>

// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/CoM.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>


namespace pg
{
class ForceContact;
class EllipseContact;

class PGData
{
public:
  struct ForceData
  {
    int bodyIndex;
    int bodyId;
    std::vector<sva::PTransformd> points;
    std::vector<sva::ForceVecd> forces;
    double mu;
  };

  struct EllipseData
  {
    int bodyIndex;  //Each ellipse is defined relatively to a Surface of a Body
    double x;     //x coord of the center
    double y;     //y coord of the center
    double theta; //Angle between the x-axis and the first axis of the ellipse
    double r1;    //First radius
    double r2;    //Second radius
    std::string print()
    {
      std::stringstream result;
      result << "ellipse = Ellipse((" << this->x << ", " << this->y << "), ";
      result << 2*this->r1 << ", " << 2*this->r2 << ", " << 180*this->theta/boost::math::constants::pi<double>() << ")\n";
      return result.str();
    }
  };

public:
  PGData(const rbd::MultiBody& mb, const Eigen::Vector3d& gravity,
         int pbSize, int qBegin, int forceBegin);

  void x(const Eigen::VectorXd& x);

  void forces(const std::vector<ForceContact>& fd);
  void ellipses(const std::vector<EllipseContact>& ed);
  void update();

  const rbd::MultiBodyConfig& mbc() const
  {
    return mbc_;
  }

  const rbd::MultiBody& multibody() const
  {
    return mb_;
  }

  const rbd::MultiBody& mb() const
  {
    return mb_;
  }

  double robotMass() const
  {
    return robotMass_;
  }

  const Eigen::Vector3d& com() const
  {
    return com_;
  }

  const Eigen::MatrixXd& comJac() const
  {
    return comJacMat_;
  }

  int pbSize() const
  {
    return pbSize_;
  }

  int qParamsBegin() const
  {
    return qBegin_;
  }

  int forceParamsBegin() const
  {
    return forceBegin_;
  }

  int ellipseParamsBegin() const
  {
    return mb_.nrParams() + nrForcePoints_*3;
  }

  int nrForcePoints() const
  {
    return nrForcePoints_;
  }

  const std::vector<ForceData>& forceDatas() const
  {
    return forceDatas_;
  }

  const std::vector<EllipseData>& ellipseDatas() const
  {
    return ellipseDatas_;
  }

  const Eigen::Vector3d& gravity() const
  {
    return gravity_;
  }

private:
  rbd::MultiBody mb_;
  rbd::MultiBodyConfig mbc_;
  double robotMass_;

  Eigen::Vector3d com_;
  rbd::CoMJacobian comJac_;
  Eigen::MatrixXd comJacMat_;

  Eigen::Vector3d gravity_;

  Eigen::VectorXd xq_, xf_;

  std::vector<ForceData> forceDatas_;
  int nrForcePoints_;

  int pbSize_;
  int qBegin_, forceBegin_;

  std::vector<EllipseData> ellipseDatas_;
};


} // namespace pg
