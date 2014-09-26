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

// associated header
#include "CoMHalfSpaceConstr.h"

// include
// std
#include <vector>

// PG
#include "PGData.h"
#include "FillSparse.h"


namespace pg
{


/*
 *                 CoMHalfSpace
 */

CoMHalfSpaceConstr::CoMHalfSpaceConstr(PGData* pgdata,
    const std::vector<Eigen::Vector3d>& O,
    const std::vector<Eigen::Vector3d>& n)
    : roboptim::DifferentiableSparseFunction(pgdata->pbSize(), O.size(), "CoMHalfSpace")
  , pgdata_(pgdata)
  , O_(O)
  , n_(n)
  , jac_(pgdata->mb())
  , eqJac_(O.size(), pgdata->mb().nrDof())
{}


CoMHalfSpaceConstr::~CoMHalfSpaceConstr()
{ }


void CoMHalfSpaceConstr::impl_compute(result_t& res, const argument_t& x) const
{
  pgdata_->x(x);
  Eigen::Vector3d C = rbd::computeCoM(pgdata_->mb(), pgdata_->mbc());

  for(std::size_t i = 0; i < O_.size(); i++)
  {
   Eigen::Vector3d OC = C - O_[i];
   res(i) = n_[i].dot(OC);
  }
}


void CoMHalfSpaceConstr::impl_jacobian(jacobian_t& jac, const argument_t& x) const
{
  pgdata_->x(x);
  jac.reserve(O_.size()*pgdata_->mb().nrDof());
  const Eigen::MatrixXd& jacMat = jac_.jacobian(pgdata_->multibody(), pgdata_->mbc());


  for(std::size_t i = 0; i < O_.size(); i++)
  {
    eqJac_.row(i) = n_[i].transpose()*jacMat;
  }
  fillSparse(eqJac_, jac, {0, pgdata_->qParamsBegin()});
}


} // pg


