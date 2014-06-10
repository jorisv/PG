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
#include "PGData.h"

// include
// std
#include <vector>

// RBDyn
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

// PG
#include "ConfigStruct.h"
#include "JacobianPatcher.h"

namespace pg
{


PGData::PGData(const rbd::MultiBody& mb, const Eigen::Vector3d& gravity,
               int pbSize, int qBegin, int forceBegin)
  : mb_(mb)
  , mbc_(mb)
  , robotMass_(0.)
  , gravity_(gravity)
  , xq_(mb.nrParams())
  , xf_()
  , nrForcePoints_(0)
  , pbSize_(pbSize)
  , qBegin_(qBegin)
  , forceBegin_(forceBegin)
  , xStamp_(0)
{
  xq_.setZero();
  mbc_.zero(mb_);
  rbd::forwardKinematics(mb_, mbc_);
  rbd::forwardVelocity(mb_, mbc_);

  for(const rbd::Body& b: mb_.bodies())
  {
    robotMass_ += b.inertia().mass();
  }
}


std::size_t PGData::x(const Eigen::VectorXd& x)
{
  if(xq_ != x.segment(qBegin_, mb_.nrParams()) ||
     xf_ != x.segment(forceBegin_, nrForcePoints_*3))
  {
    xq_ = x.segment(qBegin_, mb_.nrParams());
    xf_ = x.segment(forceBegin_, nrForcePoints_*3);
    update();
  }

   return xStamp_;
}


void PGData::forces(const std::vector<ForceContact>& forceContacts)
{
  forceDatas_.clear();
  forceDatas_.reserve(forceContacts.size());
  nrForcePoints_ = 0;
  for(const ForceContact& fc: forceContacts)
  {
    std::vector<sva::PTransformd> points(fc.points.size());
    std::vector<sva::ForceVecd> forces(fc.points.size());
    for(std::size_t i = 0; i < fc.points.size(); ++i)
    {
      points[i] = fc.points[i];
      forces[i] = sva::ForceVecd(Eigen::Vector6d::Zero());
      ++nrForcePoints_;
    }
    forceDatas_.push_back({mb_.bodyIndexById(fc.bodyId), fc.bodyId,
                           points, forces, fc.mu});
  }
  xf_.setZero(nrForcePoints_*3);
}


void PGData::ellipses(const std::vector<EllipseContact>& ellipseContacts)
{
  std::vector<PGData::EllipseData> ellipseDatas;
  ellipseDatas.reserve(ellipseContacts.size());
  for(const EllipseContact& ec: ellipseContacts)
  {
    ellipseDatas.push_back({mb_.bodyIndexById(ec.bodyId), 0., 0., 0., 0., 0.});
  }
}


void PGData::update()
{
  ++xStamp_;
  rbd::vectorToParam(xq_, mbc_.q);
  rbd::forwardKinematics(mb_, mbc_);
  patchMbc(mb_, mbc_);

  int fPos = 0;
  for(ForceData& fd: forceDatas_)
  {
    for(std::size_t i = 0; i < fd.forces.size(); ++i)
    {
      Eigen::Vector3d force(xf_[fPos + 0], xf_[fPos + 1], xf_[fPos + 2]);
      fd.forces[i] = sva::ForceVecd(Eigen::Vector3d::Zero(), force);
      fPos += 3;
    }
  }

  /*
  for(EllipseData& ed: ellipseDatas_)
  {
    construct_f()(int(x_.size()), xPos + 0, x_[xPos + 0], ed.x);
    construct_f()(int(x_.size()), xPos + 1, x_[xPos + 1], ed.y);
    construct_f()(int(x_.size()), xPos + 2, x_[xPos + 2], ed.theta);
    construct_f()(int(x_.size()), xPos + 3, x_[xPos + 3], ed.r1);
    construct_f()(int(x_.size()), xPos + 4, x_[xPos + 4], ed.r2);
    xPos += 5;
  }
  */
}


void PGData::updateKinematics(const std::vector<std::vector<double>> q)
{
  mbc_.q = q;
  rbd::forwardKinematics(mb_, mbc_);
  patchMbc(mb_, mbc_);
}


} // pg
