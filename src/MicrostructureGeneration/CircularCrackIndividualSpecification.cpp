/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CircularCrackIndividualSpecification_cpp_
#define model_CircularCrackIndividualSpecification_cpp_

#include <CircularCrackIndividualSpecification.h>

namespace model
{
    CircularCrackIndividualSpecification::CircularCrackIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("CircularCrack","Individual")
    {
        
    }

    CircularCrackIndividualSpecification::CircularCrackIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("CircularCrack","Individual",fileName)
    {
        radii=this->parser->readArray<double>("radii",true);
        if(radii.size())
        {
            nodesNumbers=this->parser->readArray<int>("nodesNumbers",true);
            meshSizes=this->parser->readArray<double>("meshSizes",true);
            centers=this->parser->readMatrix<double>("centers",radii.size(),3,true);
            normals=this->parser->readMatrix<double>("normals",radii.size(),3,true);
        }
    }
}
#endif
