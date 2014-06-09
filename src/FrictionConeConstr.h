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

// std
#include <vector>

// RBDyn
#include <RBDyn/Jacobian.h>


namespace pg
{
class PGData;

class FrictionConeConstr : public roboptim::DifferentiableSparseFunction
{
public:
  typedef typename parent_t::argument_t argument_t;

public:
  FrictionConeConstr(PGData* pgdata);
  ~FrictionConeConstr();


  void impl_compute(result_t& res, const argument_t& x) const;
  void impl_jacobian(jacobian_t& jac, const argument_t& x) const;
  void impl_gradient(gradient_t& /* gradient */,
      const argument_t& /* x */, size_type /* functionId */) const
  {
    throw std::runtime_error("NEVER GO HERE");
  }

private:
  PGData* pgdata_;
  int nrNonZero_;

  mutable std::vector<rbd::Jacobian> jacPoints_;
  mutable std::vector<Eigen::MatrixXd> jacPointsMatTmp_;
};

} // namespace pg
