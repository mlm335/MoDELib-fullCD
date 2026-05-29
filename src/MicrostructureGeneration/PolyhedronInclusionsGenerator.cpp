/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionsGenerator_cpp_
#define model_InclusionsGenerator_cpp_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PolyhedronInclusionsGenerator.h>
#include <PlaneLineIntersection.h>
#include <GmshReader.h>

namespace model
{

    bool PolyhedronInclusionsGenerator::generateSingle(MicrostructureGenerator& mg,const std::map<size_t,Eigen::Vector3d>& polyNodes,const std::map<size_t,std::vector<size_t>>& faceMap, const Eigen::Matrix<double,1,dim*dim>& eTrow, const double& vrc,const int&type)
    {
        Eigen::Matrix<double,dim,dim> eT(Eigen::Map<const Eigen::Matrix<double,dim,dim>>(eTrow.data(),dim,dim).transpose());
        mg.insertInclusion(polyNodes,faceMap, eT, vrc,type);
        return true;
    }

    PolyhedronInclusionsGenerator::PolyhedronInclusionsGenerator(const PolyhedronInclusionIndividualSpecification& spec,MicrostructureGenerator& mg) :
    /* init */ allowOverlap(false)
    {
        GmshReader mshReader(spec.mshFile);
        
//        std::cout<<spec.X0<<std::endl;
//        std::cout<<spec.F<<std::endl;
//        std::cout<<spec.eigenDistortion<<std::endl;
        
        Eigen::Vector3d nodeBary(Eigen::Vector3d::Zero());
        for(const auto& node : mshReader.nodes())
        {
            nodeBary+=node.second;
        }
        nodeBary/=mshReader.nodes().size();
        
        for(int j=0;j<spec.X0.rows();++j)
        {
            std::map<size_t,Eigen::Vector3d> scaledNodes;
            const Eigen::Vector3d x0(spec.X0.row(j));
            for(const auto& node : mshReader.nodes())
            {
                scaledNodes.emplace(node.first,spec.F*(node.second-nodeBary)+x0);
            }
            
            std::map<size_t,std::vector<size_t>> faces;
            size_t eleCounter(0);
            for(const auto& ele : mshReader.elements())
            {
                if(ele.second.type==2)
                {// 3-nodes triangle
                    faces.emplace(eleCounter,ele.second.nodeIDs);
                    eleCounter++;
                }
                
            }
            
            generateSingle(mg,scaledNodes,faces, spec.eigenDistortion, spec.velocityReductionFactor,spec.phaseID);
        }
    }
}
#endif
