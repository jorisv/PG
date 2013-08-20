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
#include <vector>

// RBDyn
#include <RBDyn/MultiBody.h>

namespace pg
{

template<typename T>
class FK
{
public:
  FK(const rbd::MultiBody& mb);

  void run(const rbd::MultiBody& mb, const std::vector<std::vector<T>>& q);

  const std::vector<sva::PTransform<T>>& bodyPosW() const
  {
    return bodyPosW_;
  }

private:
  std::vector<sva::PTransform<T>> Xt_;
  std::vector<sva::PTransform<T>> bodyPosW_;
};


// inline


template<typename T>
FK<T>::FK(const rbd::MultiBody& mb)
  : Xt_()
  , bodyPosW_(mb.nrBodies())
{
  Xt_.reserve(mb.nrJoints());
  for(const sva::PTransformd& Xt: mb.transforms())
  {
    Xt_.push_back(Xt.cast<T>());
  }
}


template<typename T>
void FK<T>::run(const rbd::MultiBody& mb, const std::vector<std::vector<T>>& q)
{
  const std::vector<rbd::Joint>& joints = mb.joints();
  const std::vector<int>& pred = mb.predecessors();
  const std::vector<int>& succ = mb.successors();

  for(std::size_t i = 0; i < joints.size(); ++i)
  {
    sva::PTransform<T> parentToSon = joints[i].pose(q[i])*Xt_[i];

    if(pred[i] != -1)
      bodyPosW_[succ[i]] = parentToSon*bodyPosW_[pred[i]];
    else
      bodyPosW_[succ[i]] = parentToSon;
  }
}


} // namespace pg
