/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CircularCrackGenerator_cpp_
#define model_CircularCrackGenerator_cpp_

#include <iostream>
#include <CircularCrackGenerator.h>

namespace model
{

CircularCrackGenerator::CircularCrackGenerator(const CircularCrackIndividualSpecification& spec,MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating individual CircularCrack"<<defaultColor<<std::endl;
    
    if(spec.radii.size()!=spec.nodesNumbers.size())
    {
        throw std::runtime_error("radii.size()="+std::to_string(spec.radii.size())+" NOT EQUAL TO nodesNumbers.size()="+std::to_string(spec.nodesNumbers.size()));
    }
    if(spec.radii.size()!=spec.meshSizes.size())
    {
        throw std::runtime_error("radii.size()="+std::to_string(spec.radii.size())+" NOT EQUAL TO meshSizes.size()="+std::to_string(spec.meshSizes.size()));
    }
    if(int(spec.radii.size())!=spec.centers.rows())
    {
        throw std::runtime_error("radii.size()="+std::to_string(spec.radii.size())+" NOT EQUAL TO centers.rows()="+std::to_string(spec.centers.rows()));
    }
    if(int(spec.radii.size())!=spec.normals.rows())
    {
        throw std::runtime_error("radii.size()="+std::to_string(spec.radii.size())+" NOT EQUAL TO normals.rows()="+std::to_string(spec.normals.rows()));
    }
    
//    double density=0.0;
    for(size_t k=0;k<spec.radii.size();++k)
    {
        try
        {
            const bool success=generateSingle(mg,spec.radii[k],spec.nodesNumbers[k],spec.meshSizes[k],spec.centers.row(k),spec.normals.row(k));
            if(success)
            {
//                density+=2.0*std::numbers::pi*spec.loopRadii[k]/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                std::cout<<"generated CircularCrack"<<std::endl;
            }
        }
        catch(const std::exception& e)
        {
            
        }
    }
}

bool CircularCrackGenerator::generateSingle(MicrostructureGenerator& mg,const double& radius,const int& nNodes,const double& meshSize, const VectorDimD& center,const VectorDimD& normal)
{
    bool success =false;
    if(radius>0)
    {
        std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(center));
        if(!found.first)
        {
            throw std::runtime_error("ShearLoopGenerator::generateSingle:: center is outside mesh.");
        }
        
        const int grainID(found.second->region->regionID);
        const auto& grain(mg.ddBase.poly.grain(grainID));
        
        Plane<3> crackPlane(center,normal);
        const VectorDimD localX1(crackPlane.L2G.col(0));
        
        std::vector<VectorDimD> crackNodePos;
        for(int k=0;k< nNodes;++k)
        {
            crackNodePos.push_back(center+Eigen::AngleAxisd(k*2.0*std::numbers::pi/nNodes,normal)*localX1*radius);
        }
        
        //    success=mg.insertJunctionLoop(loopNodePos,glidePlane,
        //                                  slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
        //                                  P0,grainID,DislocationLoopIO<3>::GLISSILELOOP);
    }
    else
    {
        std::cout<<"Skipping CircularCrack of (negative) radius "<<radius<<std::endl;
    }
    
    return success;
}

}
#endif
