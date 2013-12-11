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

#include <Eigen/Dense>

namespace pg
{

template <class problem>
struct IpoptIntermediateCallback
{
  void record(const Eigen::VectorXd & xVec, const problem &)
  {
    datas.push_back(xVec);
  }
  std::vector<Eigen::VectorXd> datas;
};

} // namespace pg
