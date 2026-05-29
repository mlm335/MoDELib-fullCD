/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationNucleation_cpp
#define model_DislocationNucleation_cpp


#include <DislocationNucleation.h>
#include <MicrostructureGenerator.h>

//#include <TypeTraits.h>
//#include <LatticeModule.h>
//#include <BCClattice.h>
//#include <FCClattice.h>
//#include <Grain.h>
//#include <TerminalColors.h>
//#include <DDtraitsIO.h>

namespace model
{

template <int dim>
DislocationNucleation<dim>::DislocationNucleation(DislocationNetwork<dim>& DN_in):
/*init*/ DN(DN_in)
{
    
}

template <int dim>
void DislocationNucleation<dim>::bulkNucleate(const int& bulkNucleationModel)
{
    switch (bulkNucleationModel)
    {
        case 1:
        {// NucleationModel #1: shear loop nucleation
            MicrostructureGenerator mg(DN.ddBase);
            for(const auto& simplex : DN.ddBase.mesh.simplices())
            {
                const auto c(simplex.second.center());
                const auto grainID(simplex.second.region->regionID);
                const auto& grain(DN.ddBase.poly.grain(grainID));
                const auto sigma(DN.microstructures.stress(c,nullptr,nullptr,nullptr));
                std::map<float,std::shared_ptr<SlipSystem>> rssMap;
                for(const auto& ss : grain->slipSystems())
                {
                    rssMap.emplace((sigma*ss.second->unitSlip).dot(ss.second->unitNormal),ss.second);
                }
                if(rssMap.begin()->first < 0.05)
                {
                    const auto& maxRssSlipSystem(rssMap.begin()->second);
                    ShearLoopIndividualSpecification spec;
                    spec.slipSystemIDs.push_back(maxRssSlipSystem->sID);
                    spec.loopRadii.push_back(50.0);
                    spec.loopCenters=c.transpose();
                    spec.loopSides.push_back(10);
                    mg.addShearLoopIndividual(spec);                    
                }
            }
            DN.addConfiguration(mg.config());
            break;
        }
            
        default:
        {
            break;
        }
    }
}

template <int dim>
void DislocationNucleation<dim>::surfaceNucleate(const int& surfaceNucleationModel)
{
    
}

template struct DislocationNucleation<3>;

}
#endif
