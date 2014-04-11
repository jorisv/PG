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
// roboptim
#include <roboptim/core/differentiable-function.hh>

// RBDyn
#include <RBDyn/Jacobian.h>


namespace pg
{
class PGData;
class BodyPositionTarget;
class BodyOrientationTarget;
class ForceContact;
class ForceContactMinimization;

class StdCostFunc : public roboptim::DifferentiableFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  StdCostFunc(PGData* pgdata, std::vector<std::vector<double>> q,
              double postureScale, double torqueScale, double forceScale, 
              double ellipseScale,
              const std::vector<BodyPositionTarget>& bodyPosTargets,
              const std::vector<BodyOrientationTarget>& bodyOriTargets,
              const std::vector<ForceContact>& forceContacts,
              const std::vector<ForceContactMinimization>& forceContactsMin);

  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_gradient(gradient_t& gradient,
      const argument_t& x, size_type /* functionId */) const throw();

private:
  struct BodyPositionTargetData
  {
    int bodyIndex;
    Eigen::Vector3d target;
    double scale;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat, jacMatFull;
  };

  struct BodyOrientationTargetData
  {
    int bodyIndex;
    Eigen::Matrix3d target;
    double scale;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat, jacMatFull;
  };

  struct ForceContactMinimizationData
  {
    std::size_t forcePos;
    std::size_t gradientPos;
    double scale;
  };

private:
  PGData* pgdata_;
  std::vector<std::vector<double>> tq_;
  double postureScale_;
  double torqueScale_;
  double forceScale_;
  double ellipseScale_;
  mutable std::vector<BodyPositionTargetData> bodyPosTargets_;
  mutable std::vector<BodyOrientationTargetData> bodyOriTargets_;
  std::vector<ForceContactMinimizationData> forceContactsMin_;
};

} // namespace pg
