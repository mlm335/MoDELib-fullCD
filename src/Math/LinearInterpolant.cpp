/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearInterpolant_cpp_
#define model_LinearInterpolant_cpp_

#include <LinearInterpolant.h>
#include <iterator>
#include <iostream>
namespace model
{

    LinearInterpolant::LinearInterpolant(const InterpolantBase::MatrixType& data,const ExtrapolationMethod& extrap):
    /* init*/ InterpolantBase(data,extrap)
    {
        
    }

LinearInterpolant::LinearInterpolant(const InterpolantBase::MapType& data,const ExtrapolationMethod& extrap):
/* init*/ InterpolantBase(data,extrap)
{
    
}

    double LinearInterpolant::_f(const double& x) const
    {
        switch (this->size())
        {
            case 0:
            {
                throw std::runtime_error("LinearInterpolant::operator() data is empty.");
                return 0.0;
                break;
            }
                
            case 1:
            {// only one point
                return this->begin()->second;
                break;
            }
                
            default:
            {// two or more points
                                
                const auto iterLower(this->lower_bound(x)); // first element in the map whose key is >= x
                const auto iterUpper(this->upper_bound(x)); // first element in the map whose key is  > x
                
                if(iterUpper==this->begin())
                {// x to the left of the domain
                    switch (this->extrapolation.type)
                    {
                        case 1:
                        {// Constant Extrapolation
                            return this->begin()->second;
                            break;
                        }
                            
                        case 2:
                        {// Periodic Extrapolation
                            InterpolantBase::MapType exMap(this->map());
                            exMap.emplace(this->map().rbegin()->first-this->extrapolation.period,this->map().rbegin()->second);
                            exMap.emplace(this->map(). begin()->first+this->extrapolation.period,this->map(). begin()->second);
                            return linearSegment(exMap.begin(),std::next(exMap.begin()),x);
                            break;
                        }
                            
                        default:
                        {// Linear Extrapolation
                            return linearSegment(this->begin(),std::next(this->begin()),x);
                            break;
                        }
                    }
                }
                else if(iterUpper==this->end())
                {// x on the right of domain
                    switch (this->extrapolation.type)
                    {
                        case 1:
                        {// Constant Extrapolation
                            return std::prev(this->end())->second;
                            break;
                        }
                            
                        case 2:
                        {// Periodic Extrapolation
                            InterpolantBase::MapType exMap(this->map());
                            exMap.emplace(this->map().rbegin()->first-this->extrapolation.period,this->map().rbegin()->second);
                            exMap.emplace(this->map(). begin()->first+this->extrapolation.period,this->map(). begin()->second);
                            return linearSegment(std::prev(exMap.end(),2),std::prev(exMap.end(),1),x);
                            break;
                        }

                        default:
                        {// Linear Extrapolation
                            return linearSegment(std::prev(this->end(),2),std::prev(this->end(),1),x);
                            break;
                        }
                    }
                }
                else
                {
                    if(iterLower==iterUpper) // x is in between breakpoints.
                    {// Because of assertion x<=xMax, iter1 cannot be dataMap.end()
                        // Moreover, iterUpper can never be dataMap.begin()
                        // Terefore it is safe to use element before iterLower
                        return linearSegment(std::prev(iterLower),iterLower,x);
                    }
                    else
                    {// x is on breakpoint iterLower. iterUpper may be dataMap.end().
                        return iterLower->second;
                    }
                }
                break;
            }
        }
    }

    double LinearInterpolant::linearSegment(const InterpolantBase::MapType::const_iterator& iter1,const InterpolantBase::MapType::const_iterator& iter2,const double& x)
    {
        return iter1->second+(iter2->second-iter1->second)/(iter2->first-iter1->first)*(x-iter1->first);
    }

} // close namespace model
#endif
