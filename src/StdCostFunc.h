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
class RobotConfig;
class RunConfig;

class StdCostFunc : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  StdCostFunc(std::vector<PGData>& pgdatas, const std::vector<RobotConfig>& robotConfigs,
              const std::vector<RunConfig>& runConfigs);

  void impl_compute(result_t& res, const argument_t& x) const;
  void impl_gradient(gradient_t& gradient,
      const argument_t& x, size_type /* functionId */) const;

private:
  struct BodyPositionTargetData
  {
    int bodyIndex;
    Eigen::Vector3d target;
    double scale;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull;
  };

  struct BodyOrientationTargetData
  {
    int bodyIndex;
    Eigen::Matrix3d target;
    double scale;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull;
  };

  struct ForceContactMinimizationData
  {
    std::size_t forcePos;
    std::size_t gradientPos;
    double scale;
  };

  struct TorqueContactMinimizationData
  {
    int bodyIndex;
    std::vector<sva::PTransformd> points;
    std::vector<Eigen::Vector3d> levers;
    Eigen::Vector3d axis;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
    Eigen::MatrixXd jacMatTmp;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull;
    std::size_t forcePos;
    std::size_t gradientPos;
    double scale;
  };

  struct NormalForceTargetData
  {
    std::size_t forcePos;
    std::size_t gradientPos;
    double target;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull;
    double scale;
  };

  struct TangentialForceMinimizationData
  {
    std::size_t forcePos;
    std::size_t gradientPos;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> jacMatFull;
    double scale;
  };

  struct RobotData
  {
    PGData* pgdata;
    std::vector<std::vector<double>> tq;
    double postureScale;
    double torqueScale;
    double forceScale;
    double ellipseScale;
    std::vector<BodyPositionTargetData> bodyPosTargets;
    std::vector<BodyOrientationTargetData> bodyOriTargets;
    std::vector<ForceContactMinimizationData> forceContactsMin;
    std::vector<TorqueContactMinimizationData> torqueContactsMin;
    std::vector<NormalForceTargetData> normalForceTargets;
    std::vector<TangentialForceMinimizationData> tanForceMin;
  };

private:
  mutable std::vector<RobotData> robotDatas_;
  double scale_;
};

} // namespace pg
