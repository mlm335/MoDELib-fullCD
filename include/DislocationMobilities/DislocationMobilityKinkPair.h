/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityKinkPair_h_
#define _model_DislocationMobilityKinkPair_h_

#include <string>
#include <PolycrystallineMaterialBase.h>
#include <DislocationMobilityBase.h>

namespace model
{
    
    struct DislocationMobilityKinkPair : public DislocationMobilityBase
    {
        static constexpr double kB_eV=8.617e-5;       // Boltzmann constant in [eV/K]

        const double kB;
        const double h;
        const double w;
        const double B0;
        const double B1;
        const double Bk;
        const double dH0;
        const double p;
        const double q;
        const double T0;
        const double tauC;
        const double a0;
        const double a1;
        const double a2;
        const double a3;
        const double conjugateAngle;

        static double sigmoid(const double & x);

        
        DislocationMobilityKinkPair(const PolycrystallineMaterialBase& material,
                                    const double& h,
                                    const double& w,
                                    const double& B0_SI,
                                    const double& B1_SI,
                                    const double& Bk_SI,
                                    const double& dH0_SI,
                                    const double& p_in,
                                    const double& q_in,
                                    const double& T0_SI,
                                    const double& tauC_SI,
                                    const double& a0_in,
                                    const double& a1_in,
                                    const double& a2_in,
                                    const double& a3_in,
                                    const double& conjugateAngle_deg) ;
        
        
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


