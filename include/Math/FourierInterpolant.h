/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FourierInterpolant_H_
#define model_FourierInterpolant_H_

#include <Eigen/Dense>

namespace model
{

    struct FourierInterpolant
    {

        const Eigen::Matrix<double,Eigen::Dynamic,2> XiYi;
        const Eigen::VectorXd W;
        
        FourierInterpolant(const Eigen::Matrix<double,Eigen::Dynamic,2>& XiFi,
                           const Eigen::Matrix<double,Eigen::Dynamic,2>& XidFi);
        
//        double p(const int& i,const double& x) const;
        double f(const double& x) const;
        
        static Eigen::VectorXd getW(const Eigen::VectorXd& Xi);
        
        
	};
    
}
#endif
