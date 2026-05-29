/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineSegmentIntersection_CPP_
#define model_LineSegmentIntersection_CPP_

#include <cfloat>
#include <tuple>
#include <map>
#include <stdexcept>
#include <Eigen/Dense>

#include <LineSegmentIntersection.h>
#include <LineLineIntersection.h>

namespace model
{
    /**********************************************************************/
    template <int dim>
    typename LineSegmentIntersection<dim>::SolutionType LineSegmentIntersection<dim>::findIntersections(const VectorDimD& A0,
                                          const VectorDimD& D0,
                                          const VectorDimD& P0,
                                          const VectorDimD& P1)
    {/*!\param[in] A0 start point of first line
      *\param[in] D0  line direction of first line
      *\param[in] P0 start point of segment
      *\param[in] P1  end point of segment
      *\returns a tuple, where the first element is an elmenet of IntersectionType (the type of intersection).
      * Second and third elements depend of the type of intersection:
      * - For COINCIDENT lines
      * - the second element is A0
      * - the third elements is A1
      * - For SKEW lines
      * - the second element is the closest point on line 1
      * - the third elements is the closest point on line 2
      * - For INCIDENT lines
      * - the second element is the point of intersection
      * - the third elements is the point of intersection
      * - For PARALLEL lines
      * - the second element is A0
      * - the third elements is A1
      * - For OUTSIDE lines
      * - the second element is the point of intersection
      * - the third elements is the point of intersection
      */
        
        
        LineLineIntersection<dim> lli(A0,D0,P0,P1-P0);
     
        switch (lli.type)
        {
            case LineLineIntersection<dim>::PARALLEL:
            {
                return SolutionType(PARALLEL,A0,0.5*(P0+P1),0.0,0.5); // TODO COMPUTE CLOSEST POINTS
                break;
            }
                
            case LineLineIntersection<dim>::COINCIDENT:
            {
                return SolutionType(COINCIDENT,A0,0.5*(P0+P1),0.0,0.5); // TODO COMPUTE CLOSEST POINTS
                break;
            }
                
            case LineLineIntersection<dim>::INCIDENT:
            {
                if(lli.u1>=0.0 && lli.u1<=1.0)
                {
                    return SolutionType(INCIDENT,lli.x0,lli.x1,lli.u0,lli.u1);
                }
                else if(lli.u1<0.0)
                {
                    return SolutionType(OUTSIDE,lli.x0,P0,lli.u0,0.0);
                }
                else //lli.u1>1.0
                {
                    return SolutionType(OUTSIDE,lli.x0,P1,lli.u0,1.0);
                }
                break;
            }
                
            default: // SKEW
            {
                return SolutionType(SKEW,A0,0.5*(P0+P1),0.0,0.5); // TODO COMPUTE CLOSEST POINTS
                break;
            }
        }
        
    }
    
    /**********************************************************************/
    template<int dim>
    LineSegmentIntersection<dim>::LineSegmentIntersection(const VectorDimD& A0,
                         const VectorDimD& D0,
                         const VectorDimD& P0,
                         const VectorDimD& P1) :
    /* init */ sol(findIntersections(A0,D0,P0,P1))
    /* init */,type(std::get<0>(sol))
    /* init */,x0(std::get<1>(sol))
    /* init */,x1(std::get<2>(sol))
    /* init */,u0(std::get<3>(sol))
    /* init */,u1(std::get<4>(sol))
    {
        
    }

    template class LineSegmentIntersection<3>;
    template class LineSegmentIntersection<2>;

} /* namespace model */
#endif
