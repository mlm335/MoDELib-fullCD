/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ShearLoopDensitySpecification_cpp_
#define model_ShearLoopDensitySpecification_cpp_

#include <ShearLoopDensitySpecification.h>

namespace model
{
    ShearLoopDensitySpecification::ShearLoopDensitySpecification():
    /* init */ MicrostructureSpecificationBase("ShearLoop","Density")
    /* init */,targetDensity_SI(0.0)
    /* init */,numberOfSides(0)
    /* init */,radiusDistributionMean_SI(0.0)
    /* init */,radiusDistributionStd_SI(0.0)
    /* init */,allowedGrainIDs({-1})
    /* init */,allowedSlipSystemIDs({-1})
    {
        
    }

    ShearLoopDensitySpecification::ShearLoopDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("ShearLoop","Density",fileName)
    /* init */,targetDensity_SI(this->parser->readScalar<double>("targetDensity_SI",true))
    /* init */,numberOfSides(targetDensity_SI>0.0? this->parser->readScalar<int>("numberOfSides",true) : 0)
    /* init */,radiusDistributionMean_SI(targetDensity_SI>0.0? this->parser->readScalar<double>("radiusDistributionMean_SI",true) : 0.0)
    /* init */,radiusDistributionStd_SI(targetDensity_SI>0.0? this->parser->readScalar<double>("radiusDistributionStd_SI",true) : 0.0)
    /* init */,allowedGrainIDs(this->parser->readArray<int>("allowedGrainIDs",true))
    /* init */,allowedSlipSystemIDs(this->parser->readArray<int>("allowedSlipSystemIDs",true))
    {
        
    }
}
#endif
