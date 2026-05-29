/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityViscousDrag_h_
#define _model_DislocationMobilityViscousDrag_h_

#include <string>
#include <PolycrystallineMaterialBase.h>
#include <DislocationMobilityBase.h>

namespace model
{
    
    struct DislocationMobilityViscousDrag : public DislocationMobilityBase
    {

        const double kB;
        const double B0;
        const double B1;
        
        DislocationMobilityViscousDrag(const PolycrystallineMaterialBase& material,
                                       const double& B0_SI,
                                       const double& B1_SI);
        
//        DislocationMobilityViscousDrag(const std::string& tag,
//                            const PolycrystallineMaterialBase& material);
        
//        std::pair<double,double> velocity(const double& tau,
//                        const double& b,
//                        const double& T) const;
//
//        std::pair<double,double> velocity(const double& tauXb,
//                        const double& T) const;
        
        double velocity(const MatrixDim& S,
                                const VectorDim& b,
                                const VectorDim& n,
                                const double& T,
                                const double& dL,
                                const double& dt,
                                const std::shared_ptr<StochasticForceGenerator>& sfg) const override ;

    };
}

#endif


