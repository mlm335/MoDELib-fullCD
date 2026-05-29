/* This file is part of PIL, the Particle Interaction Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeDirection_h_
#define model_ReciprocalLatticeDirection_h_

#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeGCD.h>

namespace model
{
    template <int dim>
    struct ReciprocalLatticeDirection :
    /* inherits */ public ReciprocalLatticeVector<dim>
    {
        typedef LatticeGCD<dim> LatticeGCDType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        ReciprocalLatticeDirection(const ReciprocalLatticeVectorType& v) ;
        ReciprocalLatticeDirection(const VectorDimI& v,const Lattice<dim>& lat);
        ReciprocalLatticeDirection(const LatticeVectorType& v1,const LatticeVectorType& v2);        
    };
    
} // end namespace
#endif
