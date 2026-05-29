/* This file is part of PIL, the Particle Interaction Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeDirection_cpp_
#define model_ReciprocalLatticeDirection_cpp_

#include <LatticeModule.h>

namespace model
{

    template <int dim>
    ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const ReciprocalLatticeVectorType& r) :
    /* delegate */ ReciprocalLatticeDirection<dim>(r.base(),r.lattice)
    {
    }

    template <int dim>
    ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const VectorDimI& v,const Lattice<dim>& lat) :
    /* base init */ ReciprocalLatticeVectorType(((v.squaredNorm()==0)? v : (v/LatticeGCD<dim>::gcd(v)).eval()),lat)
    {
    }

    template <int dim>
    ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const LatticeVectorType& v1,const LatticeVectorType& v2) :
    /* delegating */ ReciprocalLatticeDirection(ReciprocalLatticeVectorType(v1.cross(v2)))
    {
        //            assert(this->squaredNorm() && "ReciprocalLatticeDirection has Zero Norm");
    }

    template struct ReciprocalLatticeDirection<3>;
} // end namespace
#endif
