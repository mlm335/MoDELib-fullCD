/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Grain_H_
#define model_Grain_H_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Core>
#include <SimplicialMesh.h>
#include <MeshRegionObserver.h>
#include <MeshRegion.h>
#include <LatticeModule.h>
#include <GlidePlaneBase.h>
#include <PolycrystallineMaterialBase.h>
#include <SlipSystem.h>
#include <SingleCrystalBase.h>

//#include <BestRationalApproximation.h>

namespace model
{
    template <int dim>
    class GrainBoundary;
    
    template <int dim>
    class Grain : public SingleCrystalBase<dim>,
    /* base    */ public std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<dim>>>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef MeshRegion<dim> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<dim>>> GrainBoundaryContainerType;
        
    public:
        
        static constexpr double roundTol=FLT_EPSILON;
        
        const MeshRegionType& region;
        const size_t& grainID;
        
        Grain(const MeshRegionType& region_in,
              const PolycrystallineMaterialBase& material,
              const std::string& polyFile
              );
        
        const GrainBoundaryContainerType& grainBoundaries() const;
        GrainBoundaryContainerType& grainBoundaries();
        std::deque<const GlidePlaneBase*> conjugatePlaneNormal(const LatticeVectorType& B,
                                                                 const ReciprocalLatticeDirectionType& N) const;
    };
    
}
#endif
