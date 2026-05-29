/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalBase_cpp_
#define model_SingleCrystalBase_cpp_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <iterator>

#include <SingleCrystalBase.h>

namespace model
{

template<int dim>
typename SingleCrystalSelector<dim>::MatrixDim SingleCrystalSelector<dim>::latticeBasis(const PolycrystallineMaterialBase& mat)
{
    if(mat.crystalStructure=="BCC")
    {
        return BCClattice<dim>::getLatticeBasis();
    }
    else if(mat.crystalStructure=="FCC")
    {
        return FCClattice<dim>::getLatticeBasis();
    }
    else if(mat.crystalStructure=="HEX")
    {
        return HEXlattice<dim>::getLatticeBasis();
    }
    else if(mat.crystalStructure=="CubicFluorite")
    {
        return CubicFluoriteCrystal<dim>::getLatticeBasis();
    }
    else
    {
        throw std::runtime_error("SingleCrystalSelector::latticeBasis: unknown crystal structure "+mat.crystalStructure);
        return MatrixDim::Identity();
    }
}

template<int dim>
typename SingleCrystalSelector<dim>::PlaneNormalContainerType SingleCrystalSelector<dim>::planeNormals(const PolycrystallineMaterialBase& mat, const LatticeType& lat)
{
    if(mat.crystalStructure=="BCC")
    {
        return BCClattice<dim>::planeNormals(mat,lat);
    }
    else if(mat.crystalStructure=="FCC")
    {
        return FCClattice<dim>::planeNormals(mat,lat);
    }
    else if(mat.crystalStructure=="HEX")
    {
        return HEXlattice<dim>::planeNormals(mat,lat);
    }
    else if(mat.crystalStructure=="CubicFluorite")
    {
        return CubicFluoriteCrystal<dim>::planeNormals(mat,lat);
    }
    else
    {
        throw std::runtime_error("SingleCrystalSelector::planeNormals: unknown crystal structure "+mat.crystalStructure);
        return PlaneNormalContainerType();
    }
}

template<int dim>
typename SingleCrystalSelector<dim>::SlipSystemContainerType SingleCrystalSelector<dim>::slipSystems(const PolycrystallineMaterialBase& mat, const LatticeType& lat,const PlaneNormalContainerType& pnc)
{
    if(mat.crystalStructure=="BCC")
    {
        return BCClattice<dim>::slipSystems(mat,lat,pnc);
    }
    else if(mat.crystalStructure=="FCC")
    {
        return FCClattice<dim>::slipSystems(mat,lat,pnc);
    }
    else if(mat.crystalStructure=="HEX")
    {
        return HEXlattice<dim>::slipSystems(mat,lat,pnc);
    }
    else if(mat.crystalStructure=="CubicFluorite")
    {
        return CubicFluoriteCrystal<dim>::slipSystems(mat,lat,pnc);
    }
    else
    {
        throw std::runtime_error("SingleCrystalSelector::slipSystems: unknown crystal structure "+mat.crystalStructure);
        return SlipSystemContainerType();
    }
}

template<int dim>
typename SingleCrystalSelector<dim>::SecondPhaseContainerType SingleCrystalSelector<dim>::secondPhases(const PolycrystallineMaterialBase& mat, const LatticeType& lat,const PlaneNormalContainerType& pnc)
{
    if(mat.crystalStructure=="BCC")
    {
        return BCClattice<dim>::secondPhases(mat,lat,pnc);
    }
    else if(mat.crystalStructure=="FCC")
    {
        return FCClattice<dim>::secondPhases(mat,lat,pnc);
    }
    else if(mat.crystalStructure=="HEX")
    {
        return HEXlattice<dim>::secondPhases(mat,lat,pnc);
    }
    else if(mat.crystalStructure=="CubicFluorite")
    {
        return CubicFluoriteCrystal<dim>::secondPhases(mat,lat,pnc);
    }
    else
    {
        throw std::runtime_error("SingleCrystalSelector::secondPhases: unknown crystal structure "+mat.crystalStructure);
        return SecondPhaseContainerType();
    }
}

template<int dim>
SingleCrystalBase<dim>::SingleCrystalBase(const PolycrystallineMaterialBase& material,
                                          const MatrixDim& Q) :
/* init */ LatticeType(SingleCrystalSelector<dim>::latticeBasis(material),Q)
/* init */,PlaneNormalContainerType(SingleCrystalSelector<dim>::planeNormals(material,*this))
/* init */,SlipSystemContainerType(SingleCrystalSelector<dim>::slipSystems(material,*this,planeNormals()))
/* init */,SecondPhaseContainerType(SingleCrystalSelector<dim>::secondPhases(material,*this,planeNormals()))
{
}

template<int dim>
const typename SingleCrystalBase<dim>::PlaneNormalContainerType& SingleCrystalBase<dim>::planeNormals() const
{
    return *this;
}

template<int dim>
const typename SingleCrystalBase<dim>::SlipSystemContainerType& SingleCrystalBase<dim>::slipSystems() const
{
    return *this;
}

template<int dim>
const typename SingleCrystalBase<dim>::SecondPhaseContainerType& SingleCrystalBase<dim>::secondPhases() const
{
    return *this;
}

template<int dim>
const GlidePlaneBase& SingleCrystalBase<dim>:: planeBase(const size_t& planeID) const
{
    return *planeNormals().at(planeID);
}

template<int dim>
const SlipSystem& SingleCrystalBase<dim>:: slipSystem(const size_t& ssID) const
{
    return *slipSystems().at(ssID);
}

template<int dim>
const SecondPhase<dim>& SingleCrystalBase<dim>:: secondPhase(const size_t& spID) const
{
    return *secondPhases().at(spID);
}

template struct SingleCrystalBase<3>;
template struct SingleCrystalSelector<3>;

}
#endif
