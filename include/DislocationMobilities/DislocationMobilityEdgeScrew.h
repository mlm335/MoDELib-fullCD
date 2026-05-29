/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityEdgeScrew_h_
#define _model_DislocationMobilityEdgeScrew_h_

#include <DislocationMobilityBase.h>
#include <DislocationMobility.h>

namespace model
{
    struct DislocationMobilityEdgeScrew : public DislocationMobility
    {
            
        const std::shared_ptr<DislocationMobilityBase>  edgeMobility;
        const std::shared_ptr<DislocationMobilityBase>  screwMobility;


        
        DislocationMobilityEdgeScrew(const PolycrystallineMaterialBase& material,
                                     const std::shared_ptr<DislocationMobilityBase>&  edgeMobility_in,
                                     const std::shared_ptr<DislocationMobilityBase>&  screwMobility_in);
                
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg) const override;
        
//        double velocity(const MatrixDim& S,
//                        const VectorDim& b,
//                        const VectorDim& xi,
//                        const VectorDim& n,
//                        const double& T) ;
        
    };
    
}
#endif
