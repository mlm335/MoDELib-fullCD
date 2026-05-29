/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CrackMeshGenerator_H_
#define model_CrackMeshGenerator_H_


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
//#include <ShearLoopDensitySpecification.h>
#include <CircularCrackIndividualSpecification.h>


namespace model
{

    class CrackMeshGenerator
    {
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        typedef Eigen::Matrix<int,3,1> VectorDimI;

        static bool generateSingle(MicrostructureGenerator& mg,const Eigen::Matrix<double,Eigen::Dynamic,3>& nodes,const Eigen::Matrix<double,Eigen::Dynamic,3>& burgers,const Eigen::Matrix<int,Eigen::Dynamic,3>& triangles);
        
    public:
        
        CrackMeshGenerator(const CrackMeshIndividualSpecification&,MicrostructureGenerator& mg);

    };

}
#endif
