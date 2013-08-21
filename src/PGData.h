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
// std
#include <cassert>
#include <vector>

// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/MultiBody.h>

// PG
#include "FK.h"

namespace pg
{

template<typename Type>
class PGData
{
public:
  typedef typename Type::scalar_t scalar_t;
  typedef typename Type::construct_f construct_f;

public:
  PGData(const rbd::MultiBody& mb);

  void x(const Eigen::VectorXd& x);


  const FK<scalar_t>& fk() const
  {
    return fk_;
  }

  const rbd::MultiBody& multibody() const
  {
    return mb_;
  }

  int pbSize() const
  {
    return int(x_.size());
  }

private:
  rbd::MultiBody mb_;

  Eigen::VectorXd x_;
  std::vector<std::vector<scalar_t>> q_;

  FK<scalar_t> fk_;
};


// inline


template<typename Type>
PGData<Type>::PGData(const rbd::MultiBody& mb)
  : mb_(mb)
  , x_(mb.nrParams())
  , q_(mb.nrJoints())
  , fk_(mb)
{
  for(int i = 0; i < mb.nrJoints(); ++i)
  {
    q_[i].resize(mb.joint(i).params());
  }
}


template<typename Type>
void PGData<Type>::x(const Eigen::VectorXd& x)
{
  assert(x.size() == x_.size());

  if(x != x_)
  {
    x_ = x;
    int xPos = 0;
    for(int i = 0; i < mb_.nrJoints(); ++i)
    {
      for(int j = 0; j < mb_.joint(i).params(); ++j)
      {
        construct_f()(int(x_.size()), xPos, x_[xPos], q_[i][j]);
        ++xPos;
      }
    }
    fk_.run(mb_, q_);
  }
}

} // namespace pg
