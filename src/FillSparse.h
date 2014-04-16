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

// includes
// std
#include <vector>

// Eigen
#include <Eigen/Sparse>


// forward declarations
namespace rbd
{
class MultiBody;
class Jacobian;
}

namespace pg
{

typedef std::pair<int, int> jac_offset_t;

void fullJacobianSparse(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
  const Eigen::Ref<const Eigen::MatrixXd>& jacMat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset=jac_offset_t(0, 0));


void updateFullJacobianSparse(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
  const Eigen::Ref<const Eigen::MatrixXd>& jacMat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset=jac_offset_t(0, 0));


template <typename Derived>
void fillSparse(const Eigen::MatrixBase<Derived>& mat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset=jac_offset_t(0, 0))
{
  for(int row = 0; row < mat.rows(); ++row)
  {
    for(int col = 0; col < mat.cols(); ++col)
    {
      res.insert(row + offset.first, col + offset.second) = mat(row, col);
    }
  }
}


} // pg
