/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalBase_H_
#define model_SingleCrystalBase_H_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include <TerminalColors.h>
#include <LatticeModule.h>
#include <GlidePlaneBase.h>
#include <SlipSystem.h>
#include <SecondPhase.h>
#include <TextFileParser.h>
#include <SingleCrystalTraits.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <HEXlattice.h>
#include <CubicFluoriteCrystal.h>

namespace model
{

    template<int dim>
    struct SingleCrystalSelector
    {
        typedef typename TypeTraits<SingleCrystalBase<dim>>::MatrixDim MatrixDim;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::LatticeType LatticeType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::PlaneNormalContainerType PlaneNormalContainerType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::SlipSystemContainerType SlipSystemContainerType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::SecondPhaseContainerType SecondPhaseContainerType;
        
        static MatrixDim latticeBasis(const PolycrystallineMaterialBase& material);
        static PlaneNormalContainerType planeNormals(const PolycrystallineMaterialBase& material, const LatticeType& lat);
        static SlipSystemContainerType slipSystems(const PolycrystallineMaterialBase& material, const LatticeType& lat,const PlaneNormalContainerType&);
        static SecondPhaseContainerType secondPhases(const PolycrystallineMaterialBase& material, const LatticeType& lat,const PlaneNormalContainerType&);
    };

    template<int dim>
    struct SingleCrystalBase : public Lattice<dim>
    /*                  */,private TypeTraits<SingleCrystalBase<dim>>::PlaneNormalContainerType
    /*                  */,private TypeTraits<SingleCrystalBase<dim>>::SlipSystemContainerType
    /*                  */,private TypeTraits<SingleCrystalBase<dim>>::SecondPhaseContainerType
    {
        
        typedef typename TypeTraits<SingleCrystalBase<dim>>::LatticeType LatticeType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::MatrixDim MatrixDim;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::PlaneNormalContainerType PlaneNormalContainerType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::SlipSystemContainerType SlipSystemContainerType;
        typedef typename TypeTraits<SingleCrystalBase<dim>>::SecondPhaseContainerType SecondPhaseContainerType;

        SingleCrystalBase(const PolycrystallineMaterialBase& material,
                          const MatrixDim& Q);
        const PlaneNormalContainerType& planeNormals() const;
        const SlipSystemContainerType& slipSystems() const;
        const SecondPhaseContainerType& secondPhases() const;
        const GlidePlaneBase& planeBase(const size_t&) const;
        const SlipSystem& slipSystem(const size_t&) const;
        const SecondPhase<dim>& secondPhase(const size_t&) const;
    };
}
#endif
