/* This file is part of PIL, the Particle Interaction Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeDirection_h_
#define model_LatticeDirection_h_

#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeGCD.h>

namespace model
{
    template <int dim>
    struct LatticeDirection : public LatticeVector<dim>
    {
        typedef LatticeGCD<dim> LatticeGCDType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;

        LatticeDirection(const LatticeVectorType& v) ;
        LatticeDirection(const VectorDimI& v,const Lattice<dim>& lat) ;
        LatticeDirection(const ReciprocalLatticeVectorType& r1,const ReciprocalLatticeVectorType& r2);
        LatticeVectorType snapToLattice(const VectorDimD& dP) const;
        VectorDimD snapToDirection(const VectorDimD& dP) const;        
    };
    
} // end namespace
#endif
