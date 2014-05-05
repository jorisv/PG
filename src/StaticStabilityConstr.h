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
#include <roboptim/core.hh>

// RBDyn
#include <RBDyn/CoM.h>
#include <RBDyn/Jacobian.h>


namespace pg
{
class PGData;

class StaticStabilityConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  StaticStabilityConstr(PGData* pgdata);
  ~StaticStabilityConstr() throw();


  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  void computeCoM() const;

private:
  PGData* pgdata_;
  Eigen::Vector3d gravityForce_;

  mutable Eigen::Vector3d com_;
  mutable rbd::CoMJacobian comJac_;
  mutable std::size_t xStamp_;
  mutable std::vector<rbd::Jacobian> jacPoints_;
  mutable Eigen::MatrixXd jacFullMat_;
  mutable Eigen::MatrixXd T_com_fi_jac_;
  mutable Eigen::MatrixXd couple_jac_;
};

} // namespace pg
