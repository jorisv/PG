
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
#include <boost/math/constants/constants.hpp>

// PG
#include "AutoDiffFunction.h"
#include "PGData.h"

#include "ConfigStruct.h"

namespace pg
{
int sizeEllipseConstraints(std::vector<EllipseContact>& ellipseContacts)
{
  std::size_t result = 0;
  for(const EllipseContact& ec: ellipseContacts)
  {
    result += ec.targetPoints.size();
    result += ec.surfacePoints.size();
  }
  return int(result);
}

template<typename Type>
class EllipseContactConstr : public AutoDiffFunction<Type, Eigen::Dynamic>
{
public:
  typedef AutoDiffFunction<Type, Eigen::Dynamic> parent_t;
  typedef typename parent_t::scalar_t scalar_t;
  typedef typename parent_t::result_ad_t result_ad_t;
  typedef typename parent_t::argument_t argument_t;

public:
  std::string printPolygon(std::string, const std::vector<Eigen::Vector2d>&) const;
  std::string printEllipse(std::string, const typename PGData<Type>::EllipseData&) const;

  EllipseContactConstr(PGData<Type>* pgdata,
                       std::vector<EllipseContact>& ellipseContacts )
    : parent_t(pgdata, pgdata->pbSize(), 
               sizeEllipseConstraints(ellipseContacts), "EllipseContacts")
    , pgdata_(pgdata)
    , ellipseContacts_(ellipseContacts)
  {
    int contactIndex = 0;
    listTargetPointsW_.resize(ellipseContacts_.size());
    for(const EllipseContact& ec: ellipseContacts)
    {
      const std::vector<Eigen::Vector2d>& targetPoints = 
                                        ec.targetPoints;
      const sva::PTransformd& targetFrame = ec.targetFrame;
      assert(int(targetPoints.size()) > 2);

      std::vector<sva::PTransform<scalar_t> > targetPointsW(targetPoints.size());
      for(std::size_t i = 0; i < targetPoints.size(); ++i)
      {
        sva::PTransformd p(Eigen::Vector3d(targetPoints[i].x(), targetPoints[i].y(), 0.));
        targetPointsW[i] = (p*targetFrame).cast<scalar_t>();
      }

      listTargetPointsW_[contactIndex] = targetPointsW;
      ++contactIndex;
    }
  }

  ~EllipseContactConstr() throw()
  { }


  void impl_compute(result_ad_t& res, const argument_t& /* x */) const
  {
    int contactIndex = 0;
    int resIndex = 0;
    
    for(const typename PGData<Type>::EllipseData& ed: pgdata_->ellipseDatas())
    {
      int bodyIndex = pgdata_->multibody().bodyIndexById(ellipseContacts_[contactIndex].bodyId);
      const std::vector<sva::PTransform<scalar_t> >& targetPointsW = 
                                    listTargetPointsW_[contactIndex];
      const sva::PTransformd& surfaceFrame = ellipseContacts_[contactIndex].surfaceFrame;
      const std::vector<Eigen::Vector2d>& surfacePoints = 
                                    ellipseContacts_[contactIndex].surfacePoints;


      //Adding the constraint for the surfaces of the robot
      for(std::size_t i = 0; i < surfacePoints.size(); ++i)
      {
        const Eigen::Vector2d& P1 = surfacePoints[i];
        const Eigen::Vector2d& P2 = surfacePoints[(i + 1) % surfacePoints.size()];
        
        //segX and segY are the coordinates of the segment in the ellipse's frame
        scalar_t segX = P2.x()-P1.x();
        scalar_t segY = P2.y()-P1.y();
        
        //pointX and pointY are the coordinates of the vector
        //center->P1 in the ellipse's frame
        scalar_t pX = ed.x - P1.x();
        scalar_t pY = ed.y - P1.y();
        
        scalar_t c = std::cos(ed.theta);
        scalar_t s = std::sin(ed.theta);
        
        res(resIndex) = -(sqrt(pow(ed.r2*(c*segX+s*segY), 2) +
                            pow(ed.r1*(-s*segX+c*segY), 2))) + (segX*pY-segY*pX);
        ++resIndex;
      }
      //std::cout << printPolygon(std::string("robot"), surfacePoints);

      //Transformations 
      const sva::PTransform<scalar_t>& transfBtoW = pgdata_->fk().bodyPosW()[bodyIndex];
      const sva::PTransform<scalar_t>& transfStoB = surfaceFrame;
      sva::PTransform<scalar_t> transfStoW = transfStoB*transfBtoW;

      //For print purpose
      //std::vector<Eigen::Vector2d> targetPointS2d(targetPointsW.size());

      //Adding the constraint for the target surfaces of the environment
      for(std::size_t i = 0; i < targetPointsW.size(); ++i)
      {
        Eigen::Vector3<scalar_t> SurfToTargetPoint = 
                       targetPointsW[i].translation() - transfStoW.translation(); 
        scalar_t T = transfStoW.rotation().row(0).dot(SurfToTargetPoint);
        scalar_t B = transfStoW.rotation().row(1).dot(SurfToTargetPoint);
        const Eigen::Vector3<scalar_t> P1(T, B, 0);

        SurfToTargetPoint = targetPointsW[(i + 1) % 
                 targetPointsW.size()].translation() - transfStoW.translation(); 
        T = transfStoW.rotation().row(0).dot(SurfToTargetPoint);
        B = transfStoW.rotation().row(1).dot(SurfToTargetPoint);
        const Eigen::Vector3<scalar_t> P2(T, B, 0);

        //targetPointS2d[i].x() = P1.x().value();
        //targetPointS2d[i].y() = P1.y().value();

        //segX and segY are the coordinates of the segment in the ellipse's frame
        scalar_t segX = P2.x()-P1.x();
        scalar_t segY = P2.y()-P1.y();
        
        //pointX and pointY are the coordinates of the vector
        //center->P1 in the ellipse's frame
        scalar_t pX = ed.x - P1.x();
        scalar_t pY = ed.y - P1.y();
        
        scalar_t c = cos(ed.theta);
        scalar_t s = sin(ed.theta);
        
        res(resIndex) = -(sqrt(pow(ed.r2*(c*segX+s*segY), 2) +
                            pow(ed.r1*(-s*segX+c*segY), 2))) + (segX*pY-segY*pX);
        ++resIndex;
      }
      //std::cout << printPolygon(std::string("projEnv"), targetPointS2d);
      //std::cout << printEllipse(std::string("contactIndex"), ed);
      ++contactIndex;
    }    
  }


private:
  PGData<Type>* pgdata_;
  std::vector<EllipseContact>& ellipseContacts_;
  std::vector<std::vector<sva::PTransform<scalar_t> > > listTargetPointsW_;
};

template<typename Type>
std::string EllipseContactConstr<Type>::printPolygon(std::string name, const std::vector<Eigen::Vector2d>& polygon) const
{
    std::stringstream result;
    result << name << "X = [ "<< polygon[0].x();
    for(std::size_t i = 1; i < polygon.size(); ++i)
    {
      result << ", " << polygon[i].x();
    }
    result << ", " << polygon[0].x() << "]\n";
    result << name << "Y = [ " << polygon[0].y();
    for(std::size_t i = 1; i < polygon.size(); ++i)
    {
      result << ", " << polygon[i].y();
    }
    result << ", " << polygon[0].y() << "]\n";
    return result.str();
}

template<typename Type>
std::string EllipseContactConstr<Type>::printEllipse(std::string name, const typename PGData<Type>::EllipseData& ellipse) const
{
    std::stringstream result;
    result << "#Ellipse " << name << "\n";
		result << "ellipse = Ellipse((" << ellipse.x << ", " << ellipse.y << "), ";
		result << 2*ellipse.r1 << ", " << 2*ellipse.r2 << ", " << 180*ellipse.theta/boost::math::constants::pi<double>() << ")\n";
		return result.str();
}
} // namespace pg
