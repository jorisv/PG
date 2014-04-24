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
#include "FillSparse.h"

// includes
// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/Jacobian.h>


namespace pg
{


template<typename T>
void fullJacobianSparse_p(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
  const Eigen::Ref<const Eigen::MatrixXd>& jacMat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset,
  const T& set)
{
  const std::vector<int>& jointsPath = jac.jointsPath();

  for(int row = 0; row < jacMat.rows(); ++row)
  {
    int jacPos = 0;
    for(std::size_t index = 0; index < jointsPath.size(); ++index)
    {
      int i = jointsPath[index];
      int nrDof = mb.joint(i).dof();
      int posInDof = mb.jointPosInDof(i);
      for(int dof = 0; dof < nrDof; ++dof)
      {
        set(res, row + offset.first, posInDof + offset.second, jacMat(row, jacPos));
        ++posInDof;
        ++jacPos;
      }
    }
  }
}


void fullJacobianSparse(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
  const Eigen::Ref<const Eigen::MatrixXd>& jacMat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset)
{
  auto insert = [](Eigen::SparseMatrix<double, Eigen::RowMajor>& res, int row, int col, double val)
  {
    res.insert(row, col) = val;
  };
  fullJacobianSparse_p(mb, jac, jacMat, res, offset, insert);
}


void updateFullJacobianSparse(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
  const Eigen::Ref<const Eigen::MatrixXd>& jacMat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset)
{
  auto coeffRef = [](Eigen::SparseMatrix<double, Eigen::RowMajor>& res, int row, int col, double val)
  {
    res.coeffRef(row, col) = val;
  };
  fullJacobianSparse_p(mb, jac, jacMat, res, offset, coeffRef);
}


void incrementFullJacobianSparse(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
  const Eigen::Ref<const Eigen::MatrixXd>& jacMat,
  Eigen::SparseMatrix<double, Eigen::RowMajor>& res,
  const jac_offset_t& offset)
{
  auto coeffRef = [](Eigen::SparseMatrix<double, Eigen::RowMajor>& res, int row, int col, double val)
  {
    res.coeffRef(row, col) += val;
  };
  fullJacobianSparse_p(mb, jac, jacMat, res, offset, coeffRef);
}


} // pg
