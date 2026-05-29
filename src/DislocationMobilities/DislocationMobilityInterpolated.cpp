/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityInterpolated_cpp_
#define _model_DislocationMobilityInterpolated_cpp_

#include <DislocationMobilityInterpolated.h>
#include <cmath>
#include <numbers>

namespace model
{

    DislocationMobilityInterpolated::DislocationMobilityInterpolated(const PolycrystallineMaterialBase& material,
                                                                     const std::shared_ptr<InterpolantBase>& interpolant_in,
                                                                     const std::map<double,std::shared_ptr<DislocationMobilityBase>>& mobilityMap_in,
                                                                     const bool& logSpace_in):
    /* init */ DislocationMobility("DislocationMobilityInterpolated for "+material.materialName)
    /* init */,interpolant(interpolant_in)
    /* init */,mobilityMap(mobilityMap_in)
    /* init */,logSpace(logSpace_in)
    {
        if(mobilityMap.size())
        {
            if(mobilityMap.begin()->first<0.0 && mobilityMap.rbegin()->first>std::numbers::pi)
            {
                throw std::runtime_error("DislocationMobilityInterpolated: breakpoints must be between 0 and 180 deg.");
            }
        }
        else
        {
            throw std::runtime_error("DislocationMobilityInterpolated: mobilityMap is empty.");
        }
    }


    double DislocationMobilityInterpolated::velocity(const MatrixDim& S,
                                            const VectorDim& b,
                                            const VectorDim& xi,
                                            const VectorDim& n,
                                            const double& T,
                                            const double& dL,
                                            const double& dt,
                                            const std::shared_ptr<StochasticForceGenerator>& sfg) const
    {

        if(logSpace)
        {
            for(const auto& mob : mobilityMap)
            {
                const double v(mob.second->velocity(S,b,n,T,dL,dt,sfg));
                if(v<0.0)
                {
                    throw std::runtime_error("DislocationMobilityInterpolated: v<0.0 incompatible with logSpace interpolation.");
                }
                interpolant->at(mob.first)=std::log(v);
            }
        }
        else
        {
            for(const auto& mob : mobilityMap)
            {
                interpolant->at(mob.first)=mob.second->velocity(S,b,n,T,dL,dt,sfg);
            }
        }
        

        
        const VectorDim unitB(b.normalized());
        const double angleRad(std::atan2(unitB.cross(xi).dot(n), unitB.dot(xi)));
        return logSpace? std::exp(interpolant->f(angleRad)) : interpolant->f(angleRad);

    }
}
#endif
