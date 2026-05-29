/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CrackMeshIndividualSpecification_cpp_
#define model_CrackMeshIndividualSpecification_cpp_

#include <CrackMeshIndividualSpecification.h>

namespace model
{
    CrackMeshIndividualSpecification::CrackMeshIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("CircularCrack","Individual")
    {
        
    }

    CrackMeshIndividualSpecification::CrackMeshIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("CrackMesh","Individual",fileName)
    {
        nodes=this->parser->readMatrixCols<double>("nodes",3,true);
        burgers=this->parser->readMatrixCols<double>("burgers",3,true);
        triangles=this->parser->readMatrixCols<int>("triangles",3,true);
    }
}
#endif
