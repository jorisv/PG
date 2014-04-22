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

// includes
// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace pg
{

void patchMbc(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
  if(mb.joint(0).type() == rbd::Joint::Planar)
  {
    mbc.motionSubspace[0](3, 0) = -mbc.q[0][2];
    mbc.motionSubspace[0](4, 0) = mbc.q[0][1];
  }
}

} // pg
