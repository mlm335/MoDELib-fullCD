/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LagrangeInterpolant_cpp_
#define model_LagrangeInterpolant_cpp_

#include <LagrangeInterpolant.h>
#include <iterator>

namespace model
{

    LagrangeInterpolant::LagrangeInterpolant(const InterpolantBase::MatrixType& data,const ExtrapolationMethod& extrap):
    /* init*/ InterpolantBase(data,extrap)
    {
        
    }

LagrangeInterpolant::LagrangeInterpolant(const InterpolantBase::MapType& data,const ExtrapolationMethod& extrap):
/* init*/ InterpolantBase(data,extrap)
{
    
}

    double LagrangeInterpolant::_f(const double& x) const
    {
        switch (this->size())
        {
            case 0:
            {
                throw std::runtime_error("LinearInterpolant::operator() data is empty.");
                return 0.0;
                break;
            }
            
            default:
            {
                switch (this->extrapolation.type)
                {
                    case 1:
                    {// Constant Extarpolation
                        const auto iterLower(this->lower_bound(x)); // first element in the map whose key is >= x
                        const auto iterUpper(this->upper_bound(x)); // first element in the map whose key is  > x

                        if(iterUpper==this->begin())
                        {
                            return this->begin()->second;

                        }
                        else if(iterUpper==this->end())
                        {
                            return std::prev(this->end())->second;
                        }
                        else
                        {
                            return lagrangePoly(x);
                        }                        
                        break;
                    }
                        
                    case 2:
                    {// Periodic Extarpolation
                        InterpolantBase::MapType exMap(this->map());
                        exMap.emplace(this->map().rbegin()->first-this->extrapolation.period,this->map().rbegin()->second);
                        exMap.emplace(this->map(). begin()->first+this->extrapolation.period,this->map(). begin()->second);
                        return LagrangeInterpolant(exMap).f(x);
                        break;
                    }
                        
                    default:
                    {// Polynomial Extrapolation
                        return lagrangePoly(x);
                        break;
                    }
                }
                break;
            }
        }
    }

    double LagrangeInterpolant::lagrangePoly(const double& x) const
    {
        double temp(0.0);
        for (InterpolantBase::const_iterator iterI=this->begin();iterI!= this->end();++iterI)
        {
            temp+=iterI->second*p(iterI,x);
        }
        return temp;
    }

    double LagrangeInterpolant::p(const InterpolantBase::MapType::const_iterator& iterI,const double& x) const
    {
        double temp(1.0);
        for (InterpolantBase::MapType::const_iterator iterJ=this->begin();iterJ!=this->end();++iterJ)
        {
            if(iterJ!=iterI)
            {
                temp*=(x-iterJ->first)/(iterI->first-iterJ->first);
            }
        }
        return temp;
    }
}
#endif
