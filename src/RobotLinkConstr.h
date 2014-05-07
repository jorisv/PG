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
class BodyLink;

class RobotLinkConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  RobotLinkConstr(PGData* pgdata1, PGData* pgdata2,
      const std::vector<BodyLink>& linkedBodies);
  ~RobotLinkConstr() throw();

  void impl_compute(result_t& res, const argument_t& x) const throw();
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const throw();
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const throw()
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  struct LinkData
  {
    sva::PTransformd body1T, body2T;
    rbd::Jacobian jac1, jac2;
    Eigen::MatrixXd jacMat1, jacMat2;
  };

private:
  PGData* pgdata1_, *pgdata2_;
  int nrNonZero_;
  mutable std::vector<LinkData> links_;
};




} // namespace pg
