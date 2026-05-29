/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CircularCrackIndividualSpecification_H_
#define model_CircularCrackIndividualSpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct CircularCrackIndividualSpecification : public MicrostructureSpecificationBase
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        std::vector<double> radii;
        std::vector<int> nodesNumbers;
        std::vector<double> meshSizes;
        Eigen::Matrix<double,Eigen::Dynamic,3> centers;
        Eigen::Matrix<double,Eigen::Dynamic,3> normals;
        
        CircularCrackIndividualSpecification();
        CircularCrackIndividualSpecification(const std::string& fileName);
    };
}
#endif
