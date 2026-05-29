/* This file is part of PIL, the Particle Interaction Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeVector_cpp_
#define model_ReciprocalLatticeVector_cpp_

#include <iostream>
#include <LatticeModule.h>

namespace model
{

    template <int dim>
    typename ReciprocalLatticeVector<dim>::BaseType& ReciprocalLatticeVector<dim>::base()
    {
        return _base;
    }

    template <int dim>
    const typename ReciprocalLatticeVector<dim>::BaseType& ReciprocalLatticeVector<dim>::base() const
    {
        return _base;
    }

    template <int dim>
    ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const LatticeType&lat) :
/* init  */ _base(VectorDimI::Zero())
 /* init */,lattice(lat)
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    template <int dim>
    ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const VectorDimD &d,
                            /*                   */ const LatticeType&lat) :
/* init base */ _base(d2cov(d, lat))
/* init base */,lattice(lat)
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

template <int dim>
ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const VectorDimI &d,
                        /*                   */ const LatticeType&lat) :
/* init base */ _base(d)
/* init base */,lattice(lat)
{ /*!@param[in] d vector in real space
      * Constructs *this by mapping d to the lattice
      */
}

    template <int dim>
    typename ReciprocalLatticeVector<dim>::ReciprocalLatticeVectorType ReciprocalLatticeVector<dim>::operator+(const ReciprocalLatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return ReciprocalLatticeVectorType((base() + other.base()).eval(), lattice);
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::ReciprocalLatticeVectorType& ReciprocalLatticeVector<dim>::operator+=(const ReciprocalLatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::ReciprocalLatticeVectorType ReciprocalLatticeVector<dim>::operator-(const ReciprocalLatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return ReciprocalLatticeVectorType((base() - other.base()).eval(), lattice);
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::ReciprocalLatticeVectorType& ReciprocalLatticeVector<dim>::operator-=(const ReciprocalLatticeVectorType &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::ReciprocalLatticeVectorType ReciprocalLatticeVector<dim>::operator*(const long int &scalar) const
    {
        return ReciprocalLatticeVectorType((base() * scalar).eval(), lattice);
    }

    template <int dim>
    long int ReciprocalLatticeVector<dim>::dot(const LatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return base().dot(other.base());
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::LatticeDirectionType ReciprocalLatticeVector<dim>::cross(const ReciprocalLatticeVectorType &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeDirectionType(LatticeVectorType(base().cross(other.base()), lattice));
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::VectorDimD ReciprocalLatticeVector<dim>::cartesian() const
    {
        return lattice.reciprocalBasis * base().template cast<double>();
    }

    template <int dim>
    double ReciprocalLatticeVector<dim>::planeSpacing() const
    {
        return 1.0 / cartesian().norm();
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::VectorDimD ReciprocalLatticeVector<dim>::interplaneVector() const
    {
        const VectorDimD c(cartesian());
        return c / c.squaredNorm();
    }

    template <int dim>
    long int ReciprocalLatticeVector<dim>::closestPlaneIndexOfPoint(const VectorDimD &P) const
    {
        assert(base().squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        const double hd(cartesian().dot(P));
        return std::lround(hd);
    }

    template <int dim>
    long int ReciprocalLatticeVector<dim>::planeIndexOfPoint(const VectorDimD &P) const
    {
        assert(base().squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        const double hd(cartesian().dot(P));
        const long int h(std::lround(hd));
        if (fabs(hd - h) > FLT_EPSILON)
        {
            std::cout << "P=" << P.transpose() << std::endl;
            std::cout << "r=" << this->cartesian().transpose() << std::endl;
            std::cout << "hd=" << std::setprecision(15) << std::scientific << hd << std::endl;
            std::cout << "h=" << h << std::endl;
            throw std::runtime_error("P in not on a lattice plane.");
        }
        return h;
    }

    template <int dim>
    long int ReciprocalLatticeVector<dim>::planeIndexOfPoint(const LatticeVector<dim> &P) const
    {
        assert(base().squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        return dot(P);
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::VectorDimI ReciprocalLatticeVector<dim>::d2cov(const VectorDimD &d,
                            const LatticeType &lat)
    {
        const VectorDimD nd(lat.latticeBasis.transpose() * d);
        const VectorDimD rd(nd.array().round());
        if ((nd - rd).norm() > roundTol)
        {
            std::cout << "d2cov, nd=" << nd.transpose() << std::endl;
            std::cout << "d2cov, rd=" << rd.transpose() << std::endl;
            throw std::runtime_error("Input vector is not a reciprocal lattice vector");
        }
        return rd.template cast<long int>();
    }

    template <int dim>
    ReciprocalLatticeVector<dim> operator*(const long int &scalar, const ReciprocalLatticeVector<dim> &L)
    {
        return L * scalar;
    }

    template class ReciprocalLatticeVector<3>;
    template ReciprocalLatticeVector<3> operator*(const long int &scalar, const ReciprocalLatticeVector<3> &L);
} // end namespace
#endif
