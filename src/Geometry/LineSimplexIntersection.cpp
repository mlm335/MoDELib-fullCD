/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineSimplexIntersection_CPP_
#define model_LineSimplexIntersection_CPP_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>

#include <LineSimplexIntersection.h>
#include <LineSegmentIntersection.h>

namespace model
{
        
    /**********************************************************************/
    template <int dim>
    std::map<float,typename LineSimplexIntersection<dim>::VectorDim> LineSimplexIntersection<dim>::lineTriangleIntersection(const VectorDim& X,const VectorDim& s,
                                        const VectorDim& P1,const VectorDim& P2,const VectorDim& P3)
    {
        
        std::map<float,VectorDim> temp;
        
        LineSegmentIntersection<dim> lineP1P2(X,s,P1,P2);
        if(lineP1P2.type==LineSegmentIntersection<dim>::INCIDENT)
        {
            temp.emplace(lineP1P2.u1,lineP1P2.x1);
        }
        
        LineSegmentIntersection<dim> lineP2P3(X,s,P2,P3);
        if(lineP2P3.type==LineSegmentIntersection<dim>::INCIDENT)
        {
            temp.emplace(lineP2P3.u1,lineP2P3.x1);
        }
        
        LineSegmentIntersection<dim> lineP3P1(X,s,P3,P1);
        if(lineP3P1.type==LineSegmentIntersection<dim>::INCIDENT)
        {
            temp.emplace(lineP3P1.u1,lineP3P1.x1);
        }
        return temp;

//        int temp=0;
//        
//        const VectorDim v12(P2-P1);
//        const VectorDim v13(P3-P1);
//        
//        VectorDim n(v12.cross(v13)); // right-handed norm to triangle P1->P2->P3
//        const double nNorm(n.norm());
//        assert(nNorm>FLT_EPSILON && "n has zero norm");
//        n/=nNorm; // normalize
//        const double nDots(n.dot(s));
//        
//        if(std::fabs(nDots)>FLT_EPSILON) // s in not parallel to triangle
//        {
//            MatrixDim M(MatrixDim::Zero());
//            M.col(0)=v12;
//            M.col(1)=v13;
//            M.col(2)=-s;
//            const VectorDim a=M.inverse()*(X-P1);
//            
//            if(   a(0)>=0.0 && a(0)<=1.0      // intersection with trignale exists
//               && a(1)>=0.0 && a(0)+a(1)<=1.0 // intersection with trignale exists
//               && a(2)>=0.0) //  intersection along positive S-vector
//            {
//                if(nDots>=0.0)
//                {
//                    temp=-1;
//                }
//                else
//                {
//                    temp=+1;
//                }
//            }
//        }
//        //            else // s is parallel to triangle
//        //            {
//        //                assert(0 && "IMPLEMENT INTERSECTION AT INFINITY");
//        ////                if((P1-X).dot(n)>=0.0) // X is below the trianlge. Positive intersection at infinity
//        ////                {
//        ////                    temp=-1;
//        ////                }
//        ////                else // X is above the trianlge. Negative intersection at infinity
//        ////                {
//        ////                    temp=+1;
//        ////                }
//        //            }
//        
//        return temp;
    }
    
    template class LineSimplexIntersection<3>;

}
#endif
