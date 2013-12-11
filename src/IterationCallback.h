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

template <class problem_t, class solverState_t>
struct IterationCallback
{
  struct Data
  {
    Eigen::VectorXd x;
    double obj, constr_viol;
  };

  void operator()(const problem_t& /*problem*/, const solverState_t& state)
  {
    Data d;

    // we only store iteration data in regular mode
    // because some quantities don't have the same meaning
    // in other mode
    std::string mode("");
    try {
      mode = state.template getParameter<std::string>("ipopt.mode");
    }
    catch(std::out_of_range & e){ std::cerr << e.what() << std::endl;}

    if(mode == "RegularMode")
    {
      d.x = state.x();
      d.obj = state.cost().get();
      d.constr_viol = state.constraintViolation().get();
      datas.push_back(d);
    }
  }

  std::vector<Data> datas;
};

} // namespace pg
