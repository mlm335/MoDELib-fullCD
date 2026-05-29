/* This file is part of PIL, the Particle Interaction Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeVector_cpp_
#define model_LatticeVector_cpp_

#include<LatticeModule.h>

namespace model
{

    template <int dim>
    typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base()
    {
        return _base;
    }

    template <int dim>
    const typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base() const
    {
        return _base;
    }

    template <int dim>
    LatticeVector<dim>::LatticeVector(const LatticeType &lat) :
    /* init      */_base(VectorDimI::Zero())
    /* init      */,lattice(lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }

    template <int dim>
    LatticeVector<dim>::LatticeVector(const VectorDimD &d,
                                      const LatticeType &lat) :
    /* init */ _base(d2contra(d, lat))
    /* init      */,lattice(lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }

    template <int dim>
    LatticeVector<dim>::LatticeVector(const VectorDimI &n,
                                      const LatticeType& lat) :
    /* init  */ _base(n)
    /* init      */,lattice(lat)
    { /*!@param[in] d vector in real space
       * Constructs *this by mapping d to the lattice
       */
    }

    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType LatticeVector<dim>::operator+(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeVectorType((base() + other.base()).eval(), lattice);
    }

    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType& LatticeVector<dim>::operator+=(const LatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType LatticeVector<dim>::operator-(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeVectorType((base() - other.base()).eval(), lattice);
    }

    template <int dim>
    typename LatticeVector<dim>::LatticeVectorType& LatticeVector<dim>::operator-=(const LatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

    template <int dim>
    LatticeVector<dim> LatticeVector<dim>::operator*(const long int& scalar) const
    {
        return LatticeVector<dim>((base() * scalar).eval(), lattice);
    }

    template <int dim>
    long int LatticeVector<dim>::dot(const ReciprocalLatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return base().dot(other.base());
    }

    template <int dim>
    typename LatticeVector<dim>::ReciprocalLatticeDirectionType LatticeVector<dim>::cross(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return ReciprocalLatticeDirectionType(ReciprocalLatticeVectorType(base().cross(other.base()), lattice));
    }

    template <int dim>
    typename LatticeVector<dim>::VectorDimD LatticeVector<dim>::cartesian() const
    {
        return lattice.latticeBasis * base().template cast<double>();
    }

    template <int dim>
    typename LatticeVector<dim>::VectorDimI LatticeVector<dim>::d2contra(const VectorDimD &d,
                                                                         const LatticeType &lat)
    {
        const VectorDimD nd(lat.reciprocalBasis.transpose() * d);
        const VectorDimD rd(nd.array().round());
        if ((nd - rd).norm() > roundTol)
        {
            std::cout << "d2contra, nd=" << nd.transpose() << std::endl;
            std::cout << "d2contra, rd=" << rd.transpose() << std::endl;
            throw std::runtime_error("Input vector is not a lattice vector");
        }
        return rd.template cast<long int>();
    }

    template<int dim>
    LatticeVector<dim> operator*(const long int& scalar, const LatticeVector<dim>& L)
    {
        return L*scalar;
    }

    template class LatticeVector<3>;
    template LatticeVector<3>operator*(const long int& scalar, const LatticeVector<3>& L);

} // end namespace
#endif
