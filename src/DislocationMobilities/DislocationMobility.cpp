/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobility_cpp_
#define _model_DislocationMobility_cpp_

#include <DislocationMobility.h>

namespace model
{

    

    DislocationMobility::DislocationMobility(const std::string& name_in) :
    /* init */ name(name_in)
    {
        std::cout<<greenBoldColor<<"Creating "<<name<<defaultColor<<std::endl;
    }

    double DislocationMobility::velocity(const MatrixDim& S,
                const VectorDim& b,
                const VectorDim& xi,
                const VectorDim& n,
                const double& T) const
    {
        return velocity(S,b,xi,n,T,0.0,0.0,nullptr);
    }

}
#endif
