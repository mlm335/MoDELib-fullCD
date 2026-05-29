/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CrackMeshIndividualSpecification_H_
#define model_CrackMeshIndividualSpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct CrackMeshIndividualSpecification : public MicrostructureSpecificationBase
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        Eigen::Matrix<double,Eigen::Dynamic,3> nodes;
        Eigen::Matrix<double,Eigen::Dynamic,3> burgers;
        Eigen::Matrix<int,Eigen::Dynamic,3> triangles;
        
        CrackMeshIndividualSpecification();
        CrackMeshIndividualSpecification(const std::string& fileName);
    };
}
#endif
