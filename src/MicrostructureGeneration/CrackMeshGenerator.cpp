/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CrackMeshGenerator_cpp_
#define model_CrackMeshGenerator_cpp_

#include <iostream>
#include <CrackMeshGenerator.h>

namespace model
{

CrackMeshGenerator::CrackMeshGenerator(const CrackMeshIndividualSpecification& spec,MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating individual CrackMesh"<<defaultColor<<std::endl;
    
    if(spec.nodes.rows()!=spec.burgers.rows())
    {
        throw std::runtime_error("nodes.rows()="+std::to_string(spec.nodes.rows())+" NOT EQUAL TO burgers.rows()="+std::to_string(spec.burgers.rows()));
    }
    
    const bool success=generateSingle(mg,spec.nodes,spec.burgers,spec.triangles);
}


bool CrackMeshGenerator::generateSingle(MicrostructureGenerator& mg,const Eigen::Matrix<double,Eigen::Dynamic,3>& nodes,const Eigen::Matrix<double,Eigen::Dynamic,3>& burgers,const Eigen::Matrix<int,Eigen::Dynamic,3>& triangles)
{
    bool success =false;
    if(nodes.rows()>2 && triangles.rows()>0)
    {
        // Convert to std::vectors
        std::vector<VectorDimD> nodesV;
        std::vector<VectorDimD> burgersV;
        std::vector<VectorDimI> trianglesV;

        for(int k=0;k<nodes.rows();++k)
        {
            nodesV.emplace_back(nodes.row(k));
            burgersV.emplace_back(burgers.row(k));
        }
        
        for(int k=0;k<triangles.rows();++k)
        {
            trianglesV.emplace_back(triangles.row(k));
        }
        
        success=mg.insertCrackMesh(nodesV,burgersV,trianglesV);

        

    }

    return success;
}

}
#endif
