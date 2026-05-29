/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Grain_cpp_
#define model_Grain_cpp_

#include <Grain.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <HEXlattice.h>

namespace model
{

    template <int dim>
    Grain<dim>::Grain(const MeshRegionType& region_in,
                      const PolycrystallineMaterialBase& material,
                      const std::string& polyFile
                      ) :
    /* init */SingleCrystalBase<dim>(material,TextFileParser(polyFile).readMatrix<double>("C2G"+std::to_string(region_in.regionID),dim,dim,true))
    /* init */,region(region_in)
    /* init */,grainID(region.regionID) // remove grain ID, use lattice.sID
    {
        if(this->sID!=region.regionID)
        {
            std::cout<<"this->sID="<<this->sID<<std::endl;
            std::cout<<"region.regionID="<<region.regionID<<std::endl;
            throw std::runtime_error("grain.sID!=region.regionID.");
        }
        std::cout<<"  # planeNormals="<<this->planeNormals().size()<<std::endl;
        std::cout<<"  # slipSystems="<<this->slipSystems().size()<<std::endl;
        std::cout<<"  # secondPhases="<<this->secondPhases().size()<<std::endl;
    }

    template <int dim>
    const typename Grain<dim>::GrainBoundaryContainerType& Grain<dim>::grainBoundaries() const
    {
        return *this;
    }

    template <int dim>
    typename Grain<dim>::GrainBoundaryContainerType& Grain<dim>::grainBoundaries()
    {
        return *this;
    }

    template <int dim>
    std::deque<const GlidePlaneBase*> Grain<dim>::conjugatePlaneNormal(const LatticeVectorType& B,
                                                                         const ReciprocalLatticeDirectionType& N) const
    {
        std::deque<const GlidePlaneBase*> temp;
        if(B.dot(N)==0) // not sessile
        {
            for (const auto& planeNormal : this->planeNormals())
            {
                if(	 B.dot(*planeNormal.second)==0 && N.cross(*planeNormal.second).base().squaredNorm()>0)
                {
                    temp.push_back(planeNormal.second.get());
                }
            }
        }
        return temp;
    }

    template class Grain<3>;
}
#endif
