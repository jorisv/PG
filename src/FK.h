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

  void init(int nrVars);
  void run(const rbd::MultiBody& mb, const std::vector<std::vector<T>>& q);

  const std::vector<sva::PTransform<T>>& bodyPosW() const
  {
    return bodyPosW_;
  }

  const std::vector<sva::PTransform<T>>& parentToSon() const
  {
    return parentToSon_;
  }

private:
  std::vector<sva::PTransform<T>> Xt_;
  std::vector<sva::PTransform<T>> bodyPosW_;
  std::vector<sva::PTransform<T>> parentToSon_;
};


// inline


template<typename T>
FK<T>::FK(const rbd::MultiBody& mb)
  : Xt_()
  , bodyPosW_(mb.nrBodies())
  , parentToSon_(mb.nrJoints())
{
  Xt_.reserve(mb.nrJoints());
  for(const sva::PTransformd& Xt: mb.transforms())
  {
    Xt_.push_back(Xt.cast<T>());
  }
}


template<typename T>
void FK<T>::init(int nrVars)
{
  for(sva::PTransform<T>& Xt: Xt_)
  {
    for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
      {
        Xt.rotation()(i,j).derivatives().setZero(nrVars);
      }
      Xt.translation()(i).derivatives().setZero(nrVars);
    }
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
    parentToSon_[i] = joints[i].pose(q[i])*Xt_[i];

    if(pred[i] != -1)
      bodyPosW_[succ[i]] = parentToSon_[i]*bodyPosW_[pred[i]];
    else
      bodyPosW_[succ[i]] = parentToSon_[i];
  }
}


} // namespace pg
