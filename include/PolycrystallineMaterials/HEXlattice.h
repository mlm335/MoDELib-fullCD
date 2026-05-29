/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HEXlattice_H_
#define model_HEXlattice_H_

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include <LatticeModule.h>
#include <SlipSystem.h>
#include <PolycrystallineMaterialBase.h>
#include <RationalLatticeDirection.h>
#include <SingleCrystalTraits.h>
#include <DislocationMobilitySelector.h>

namespace model
{

    template<int dim>
    struct HEXlattice
    {
        
    };

    template<>
    struct HEXlattice<3>
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::LatticeType LatticeType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::MatrixDim MatrixDim;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::PlaneNormalContainerType PlaneNormalContainerType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::SlipSystemContainerType SlipSystemContainerType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::SecondPhaseContainerType SecondPhaseContainerType;

        static Eigen::Matrix<double,dim,dim> getLatticeBasis();
        static PlaneNormalContainerType planeNormals(const PolycrystallineMaterialBase& material,const LatticeType& lat);
        static SlipSystemContainerType slipSystems(const PolycrystallineMaterialBase& material,const LatticeType& lat,const PlaneNormalContainerType& plN);
        static SecondPhaseContainerType secondPhases(const PolycrystallineMaterialBase& material,const LatticeType& lat,const PlaneNormalContainerType& plN);
    };

}
#endif

