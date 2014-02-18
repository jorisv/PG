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

// Eigen
#include <Eigen/Core>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/MultiBody.h>

// PG
#include "FK.h"


namespace pg
{

template<typename T>
class ID
{
public:
  ID(const rbd::MultiBody& mb, const Eigen::Vector3d& gravity);

  void run(const rbd::MultiBody& mb,
           const std::vector<sva::PTransform<T>>& bodyPosW,
           const std::vector<sva::PTransform<T>>& parentToSon,
           const std::vector<sva::ForceVec<T>>& forces);

  const std::vector<sva::MotionVec<T>>& bodyAcc() const
  {
    return bodyAcc_;
  }

  const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& torque() const
  {
    return t_;
  }

private:
  sva::MotionVec<T> gravity_;
  std::vector<sva::RBInertia<T>> I_;
  std::vector<Eigen::Matrix<double, 6, Eigen::Dynamic>> S_;
  std::vector<sva::RBInertia<T>> Ic_;
  std::vector<sva::MotionVec<T>> bodyAcc_;
  std::vector<sva::ForceVec<T>> bodySupFor_;
  std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> t_;
};


// inline


template<typename T>
ID<T>::ID(const rbd::MultiBody& mb, const Eigen::Vector3d& gravity)
  : gravity_(Eigen::Vector3<T>::Zero(), gravity.cast<T>())
  , I_(mb.nrBodies())
  , S_(mb.nrJoints())
  , Ic_(mb.nrBodies())
  , bodyAcc_(mb.nrBodies())
  , bodySupFor_(mb.nrBodies())
  , t_(mb.nrJoints())
{
  for(int i = 0; i < mb.nrBodies(); ++i)
  {
    I_[i] = mb.body(i).inertia().cast<T>();
  }

  for(int i = 1; i < mb.nrJoints(); ++i)
  {
    S_[i] = mb.joint(i).motionSubspace();
    t_[i].resize(mb.joint(i).dof(), 1);
  }
}


template<typename T>
void ID<T>::run(const rbd::MultiBody& mb,
                const std::vector<sva::PTransform<T>>& /* bodyPosW */,
                const std::vector<sva::PTransform<T>>& parentToSon,
                const std::vector<sva::ForceVec<T>>& forces)
{
  const std::vector<int>& parents = mb.parents();

  bodyAcc_[0] = parentToSon[0]*gravity_;
  for(int i = 1; i < mb.nrBodies(); ++i)
  {
    bodyAcc_[i] = parentToSon[i]*bodyAcc_[parents[i]];
    Ic_[i] = I_[i];
    bodySupFor_[i] = I_[i]*bodyAcc_[i] - forces[i];
  }

  Ic_[0] = I_[0];
  bodySupFor_[0] = I_[0]*bodyAcc_[0] - forces[0];

  for(int i = mb.nrBodies() - 1; i > 0; --i)
  {
    int parent = parents[i];
    Ic_[parent] += parentToSon[i].transMul(Ic_[i]);
    bodySupFor_[parent] += parentToSon[i].transMul(bodySupFor_[i]);
  }

  bodyAcc_[0] = sva::MotionVec<T>(-(Ic_[0].matrix().ldlt().solve(bodySupFor_[0].vector())));

  for(int i = 1; i < mb.nrBodies(); ++i)
  {
    bodyAcc_[i] = parentToSon[i]*bodyAcc_[parents[i]];
    t_[i] = S_[i].transpose()*(Ic_[i]*bodyAcc_[i] + bodySupFor_[i]).vector();
  }
}


} // namespace pg
