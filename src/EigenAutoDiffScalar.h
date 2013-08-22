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

// Eigen
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>

namespace pg
{

struct eigen_ad
{
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> derivative_t;
  typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;

  struct construct_f
  {
    void operator()(int nrX, int pos, double value, scalar_t& diff)
    {
      diff = scalar_t(value, nrX, pos);
    }
  };
};

} // namespace pg


namespace std
{

pg::eigen_ad::scalar_t pow(const pg::eigen_ad::scalar_t& v, int p)
{
  return Eigen::pow(v, p);
}

pg::eigen_ad::scalar_t sqrt(const pg::eigen_ad::scalar_t& v)
{
  return Eigen::sqrt(v);
}

pg::eigen_ad::scalar_t sin(const pg::eigen_ad::scalar_t& v)
{
  return Eigen::sin(v);
}

pg::eigen_ad::scalar_t cos(const pg::eigen_ad::scalar_t& v)
{
  return Eigen::cos(v);
}

pg::eigen_ad::scalar_t acos(const pg::eigen_ad::scalar_t& v)
{
  return Eigen::acos(v);
}

pg::eigen_ad::scalar_t log(const pg::eigen_ad::scalar_t& v)
{
  return Eigen::log(v);
}

} // namespace std
