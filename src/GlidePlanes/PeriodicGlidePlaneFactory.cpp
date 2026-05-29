/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlaneFactory_cpp_
#define model_PeriodicGlidePlaneFactory_cpp_

#include <Eigen/Dense>
#include <PeriodicGlidePlaneFactory.h>

namespace model
{

    template<int dim>
    PeriodicGlidePlaneFactory<dim>::PeriodicGlidePlaneFactory(const Polycrystal<dim>& poly_in,
                              GlidePlaneFactory<dim>& glidePlaneFactory_in) :
    /* init */ poly(poly_in)
    /* init */,glidePlaneFactory(glidePlaneFactory_in)
    {
        for(const auto& region : poly.mesh.regions())
        {
            for(const auto& face : region.second->faces())
            {
                if(face.second->periodicFacePair.second)
                {
                    if(  int(region.second->regionID)==face.second->periodicFacePair.second->regionIDs.first
                       ||int(region.second->regionID)==face.second->periodicFacePair.second->regionIDs.second)
                    {// periodic translation re-enters into the same grain
                        const auto& grain(poly.grains.at(region.second->regionID));
                        std::cout<<"MeshFace "<<face.first<<" checking if periodicShift is a LatticeVector "<<std::flush;
                        const auto lv(grain->latticeVector(face.second->periodicFacePair.first));
                        std::cout<<"LatticeVector"<<"="<<lv.base().transpose()<<std::endl;
                    }
                }
            }
        }
    }
    
    template<int dim>
    typename PeriodicGlidePlaneFactory<dim>::PeriodicGlidePlaneSharedPtrType  PeriodicGlidePlaneFactory<dim>::get(const GlidePlaneKeyType& temp)
    {
        return BaseType::getFromKey(temp);
    }
    
    template<int dim>
    typename PeriodicGlidePlaneFactory<dim>::BaseType&  PeriodicGlidePlaneFactory<dim>::periodicGlidePlanes()
    {
        return *this;
    }
    
    template<int dim>
    const typename PeriodicGlidePlaneFactory<dim>::BaseType&  PeriodicGlidePlaneFactory<dim>::periodicGlidePlanes() const
    {
        return *this;
    }
    
    template struct PeriodicGlidePlaneFactory<3>;
}
#endif
