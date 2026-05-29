/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityInterpolated_h_
#define _model_DislocationMobilityInterpolated_h_

#include <map>
#include <DislocationMobilityBase.h>
#include <DislocationMobility.h>
#include <InterpolantBase.h>

namespace model
{
    struct DislocationMobilityInterpolated : public DislocationMobility
    {
        const std::shared_ptr<InterpolantBase> interpolant;
        const std::map<double,std::shared_ptr<DislocationMobilityBase>> mobilityMap; // map of <angle, mobility>
        const bool logSpace;
                
        DislocationMobilityInterpolated(const PolycrystallineMaterialBase& material,
                                        const std::shared_ptr<InterpolantBase>& interpolant_in,
                                        const std::map<double,std::shared_ptr<DislocationMobilityBase>>& mobilityMap_in,
                                        const bool& logSpace_in=false);
                
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg) const override;
        
        
    };
    
}
#endif
