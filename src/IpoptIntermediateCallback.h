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
#include <roboptim/core/plugin/ipopt.hh>

// Ipopt
#include <IpIpoptData.hpp>
#include <IpIteratesVector.hpp>
#include <IpDenseVector.hpp>
#include <IpIpoptCalculatedQuantities.hpp>

namespace pg
{

struct IpoptIntermediateCallback : public roboptim::UserIntermediateCallback
{
  struct Data
  {
    Eigen::VectorXd x;
    double obj, dual_inf, constr_viol, complem, overallError;
  };

  bool operator()(Ipopt::AlgorithmMode mode,
    int /* iter */, double /* obj_value */,
    double /* inf_pr */, double /* inf_du */,
    double mu, double /* d_norm */,
    double /* regularization_size */,
    double /* alpha_du */, double /* alpha_pr */,
    int /* ls_trials */,
    const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq)
  {
    Data d;

    const Ipopt::DenseVector* x;

    // we only store iteration data in regular mode
    // because some quantities don't have the same meaning
    // in other mode
    if(mode == Ipopt::RegularMode)
    {
      x = dynamic_cast<const Ipopt::DenseVector*>(
        GetRawPtr(ip_data->curr()->x()));

      d.x = Eigen::Map<const Eigen::VectorXd>(x->Values(), x->Dim());
      d.obj = ip_cq->curr_f();
      d.dual_inf = ip_cq->unscaled_curr_dual_infeasibility(Ipopt::NORM_MAX);
      d.constr_viol = ip_cq->unscaled_curr_nlp_constraint_violation(Ipopt::NORM_MAX);
      d.complem = ip_cq->unscaled_curr_complementarity(mu, Ipopt::NORM_MAX);
      d.overallError = ip_cq->curr_nlp_error();
      datas.push_back(d);
    }

    return true;
  }

  std::vector<Data> datas;
};

} // namespace pg
