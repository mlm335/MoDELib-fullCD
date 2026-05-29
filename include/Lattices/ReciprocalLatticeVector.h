/* This file is part of PIL, the Particle Interaction Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ReciprocalLatticeVector_h_
#define model_ReciprocalLatticeVector_h_

#include <iostream>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <Lattice.h>

namespace model
{
    
    
    template <int dim>
    class ReciprocalLatticeVector //: public Eigen::Matrix<long int,dim,1>
    {
        static_assert(dim>0,"dim must be > 0.");

        Eigen::Matrix<long int,dim,1> _base;

    public:

        typedef Eigen::Matrix<long int,dim,1> BaseType;
        typedef Lattice<dim> LatticeType;

        BaseType& base();
        const BaseType& base() const;
                
        static constexpr double roundTol=FLT_EPSILON;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;

        const LatticeType& lattice;
        
        ReciprocalLatticeVector(const LatticeType& lat);
        ReciprocalLatticeVector(const VectorDimD& d,const LatticeType& lat) ;
        ReciprocalLatticeVector(const VectorDimI& d,const LatticeType& lat) ;
        ReciprocalLatticeVectorType operator+(const ReciprocalLatticeVectorType& other) const;
        ReciprocalLatticeVectorType& operator+=(const ReciprocalLatticeVectorType& other);
        ReciprocalLatticeVectorType operator-(const ReciprocalLatticeVectorType& other) const;
        ReciprocalLatticeVectorType& operator-=(const ReciprocalLatticeVectorType& other);
        ReciprocalLatticeVectorType operator*(const long int& scalar) const;
        long int dot(const LatticeVectorType& other) const;
        LatticeDirectionType cross(const ReciprocalLatticeVectorType& other) const;
        VectorDimD cartesian() const;
        static VectorDimI d2cov(const VectorDimD& d,const LatticeType& lat);
        double planeSpacing() const;
        VectorDimD interplaneVector() const;
        long int closestPlaneIndexOfPoint(const VectorDimD& P) const;
        long int planeIndexOfPoint(const VectorDimD& P) const ;
        long int planeIndexOfPoint(const LatticeVector<dim>& P) const;
    };
    
    template<int dim>
    ReciprocalLatticeVector<dim> operator*(const long int& scalar, const ReciprocalLatticeVector<dim>& L);
    
} // end namespace
#endif
