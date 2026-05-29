/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InterpolantBase_cpp_
#define model_InterpolantBase_cpp_

#include <iostream>
#include <algorithm>
#include <InterpolantBase.h>

namespace model
{

ExtrapolationMethod::ExtrapolationMethod():
/* init */ type(0)
/* init */,period(0.0)
{
    
}

ExtrapolationMethod::ExtrapolationMethod(const int& t,const double& p):
/* init */ type(t)
/* init */,period(std::fabs(p))
{
    
}
    
    InterpolantBase::InterpolantBase(const MapType& data,const ExtrapolationMethod& extrap) :
    /* init */ MapType(data)
    /* init */,extrapolation(extrap)
    {/*! Constructure reads data abd initializes dataMap.
      */
        if(map().size()==0)
        {
            throw std::runtime_error("InterpolantBase: sortedData is empty.");
        }
        
        if(extrapolation.type==2)
        {
            if(extrapolation.period<xMax()-xMin())
            {
                throw std::runtime_error("InterpolantBase: extrapolation.period < xMax-xMin.");
            }
        }
    }

    InterpolantBase::InterpolantBase(const InterpolantBase::MatrixType& data,const ExtrapolationMethod& extrap) :
    /* init */ InterpolantBase(matrix2Map(data),extrap)
    {
    }

    InterpolantBase::~InterpolantBase()
    {
        
    }

    typename InterpolantBase::MapType InterpolantBase::matrix2Map(const MatrixType& data)
    {
        MapType temp;
        for(int r=0;r<data.rows();++r)
        {
            const auto success(temp.emplace(data(r,0),data(r,1)));
            if(!success.second)
            {
                std::cout<<"x="<<data(r,0)<<", y="<<data(r,1)<<std::endl;
                throw std::runtime_error("InterpolantBase::convertData found non-unique x values.");
            }
        }
        return temp;
    }

typename InterpolantBase::MatrixType InterpolantBase::map2Matrix(const MapType& data)
{
    MatrixType temp(data.size(),2);
    int k(0);
    for(const auto& pair : data)
    {
        temp(k,0)=pair.first;
        temp(k,1)=pair.second;
        k++;
    }
    return temp;
}

    double InterpolantBase::xMin() const
    {
        return this->begin()->first;
    }

    double InterpolantBase::xMax() const
    {
        return this->rbegin()->first;
    }

    double InterpolantBase::f(const double& x) const
    {
        if(extrapolation.type==2)
        {
            return _f(std::fmod(x-xMin(),extrapolation.period)+xMin()+(x-xMin()<0.0? extrapolation.period : 0.0));
        }
        else
        {
            return _f(x);
        }
    }

    double InterpolantBase::operator()(const double& x) const
    {
        return f(x);
    }

//double InterpolantBase::atPeriodic(const double& x,const double& a,const double& b) const
//{
//    const double xL(std::min(a,b));
//    const double xH(std::max(a,b));
//    
//    if(xL>xMin() || xH<xMax())
//    {
//        throw std::runtime_error("InterpolantBase::atPeriodic periodic limits must include data range.");
//        return 0.0;
//    }
//        
//    const double y(std::fmod(x-xL, xH-xL)+xL); SOMETHING WRONG HERE
//    
//    return this->operator()(y);
//}



    InterpolantBase::MapType& InterpolantBase::map()
    {
        return *this;
    }

    const InterpolantBase::MapType& InterpolantBase::map() const
    {
        return *this;
    }

} // close namespace model
#endif
