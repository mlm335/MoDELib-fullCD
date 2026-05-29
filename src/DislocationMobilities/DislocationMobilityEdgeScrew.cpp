/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityEdgeScrew_cpp_
#define _model_DislocationMobilityEdgeScrew_cpp_

#include <DislocationMobilityEdgeScrew.h>

namespace model
{

DislocationMobilityEdgeScrew::DislocationMobilityEdgeScrew(const PolycrystallineMaterialBase& material,
                                                           const std::shared_ptr<DislocationMobilityBase>&  edgeMobility_in,
                             const std::shared_ptr<DislocationMobilityBase>&  screwMobility_in):
/* init */ DislocationMobility("DislocationMobilityEdgeScrew for "+material.materialName)
/* init */,edgeMobility(edgeMobility_in)
/* init */,screwMobility(screwMobility_in)
{
    
}


    double DislocationMobilityEdgeScrew::velocity(const MatrixDim& S,
                                            const VectorDim& b,
                                            const VectorDim& xi,
                                            const VectorDim& n,
                                            const double& T,
                                            const double& dL,
                                            const double& dt,
                                            const std::shared_ptr<StochasticForceGenerator>& sfg) const
    {

        const double cos2=std::pow(b.normalized().dot(xi),2);
        const double vs(screwMobility? screwMobility->velocity(S,b,n,T,dL,dt,sfg) : 0.0);
        const double ve( edgeMobility?  edgeMobility->velocity(S,b,n,T,dL,dt,sfg) : 0.0);
        return vs*cos2+ve*(1.0-cos2);
    }
}
#endif
