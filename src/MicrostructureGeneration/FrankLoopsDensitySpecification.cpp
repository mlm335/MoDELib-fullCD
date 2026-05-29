/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsDensitySpecification_cpp_
#define model_FrankLoopsDensitySpecification_cpp_

#include <FrankLoopsDensitySpecification.h>

namespace model
{
    FrankLoopsDensitySpecification::FrankLoopsDensitySpecification():
    /* init */ MicrostructureSpecificationBase("PeriodicDipole","Density")
    /* init */,targetDensity(0.0)
    /* init */,numberOfSides(0)
    /* init */,radiusDistributionMean(0.0)
    /* init */,radiusDistributionStd(0.0)
    /* init */,areVacancyLoops(true)
    /* init */,allowedGrainIDs({-1})
    /* init */,allowedPlaneIDs({-1})
    {
        
    }

    FrankLoopsDensitySpecification::FrankLoopsDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("FrankLoops","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity",true))
    /* init */,numberOfSides(targetDensity>0.0? this->parser->readScalar<int>("numberOfSides",true) : 0)
    /* init */,radiusDistributionMean(targetDensity>0.0? this->parser->readScalar<double>("radiusDistributionMean_SI",true) : 0.0)
    /* init */,radiusDistributionStd(targetDensity>0.0? this->parser->readScalar<double>("radiusDistributionStd_SI",true) : 0.0)
    /* init */,areVacancyLoops(targetDensity>0.0? this->parser->readScalar<int>("areVacancyLoops",true) : 1)
    /* init */,allowedGrainIDs(this->parser->readArray<int>("allowedGrainIDs",true))
    /* init */,allowedPlaneIDs(this->parser->readArray<int>("allowedPlaneIDs",true))
    {
        
    }
}
#endif
