/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationNETWORKTRAITS_H_
#define model_DislocationNETWORKTRAITS_H_

#include <Eigen/Dense>
#include <TypeTraits.h>

//#include <Quadrature.h>
//#include <QuadratureDynamic.h>
//#include <QuadPowDynamic.h>
#include <LatticeModule.h>
//#include <SplineBase.h>
#include <GlidePlaneModule.h>

namespace model
{
    
    
    
    /************************************************************/
    /*	Class Predeclarations ***********************************/
    template <int dim>
    class DislocationNetwork;
    
    template <int dim>
    class DislocationLoop;
    
    template <int dim>
    class DislocationLoopNode;
    
    template <int dim>
    class DislocationLoopLink;
    
    template <int dim>
    class DislocationNode;
    
    template <int dim>
    class DislocationSegment;
    
    template <typename DislocationNetworkType>
    class DislocationNetworkRemesh;
    
    template <typename DislocationNetworkType>
    class DislocationCrossSlip;
    
    template <typename DislocationNetworkType>
    class DislocationJunctionFormation;
    
    /********************************************************************/
    /*	DislocationNetworkTraitsBase: a base class for Dislocation Network Traits */
    template <int _dim>
    struct DislocationNetworkTraitsBase
    {
        static constexpr int dim=_dim;
        typedef DislocationNetwork   <dim>	LoopNetworkType;
        typedef DislocationLoop      <dim> LoopType;
        typedef DislocationLoopNode  <dim>	LoopNodeType;
        typedef DislocationLoopLink  <dim> LoopLinkType;
        typedef DislocationNode      <dim> NetworkNodeType;
        typedef DislocationSegment   <dim>	NetworkLinkType;
        
        
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef RationalLatticeDirection<dim> RationalLatticeDirectionType;
        
        typedef RationalLatticeDirectionType                         FlowType;
        //        typedef double                         FlowType;
        typedef Eigen::Matrix<double,dim,1>                         VectorDim;
        typedef Eigen::Matrix<double,dim-1,1>                         VectorLowerDim;
        typedef Eigen::Matrix<double,dim,dim>                       MatrixDim;
        typedef Grain<dim>                       GrainType;
        typedef GlidePlane<dim>                       GlidePlaneType;
        typedef PeriodicGlidePlane<dim>                       PeriodicGlidePlaneType;
        
        //        static constexpr FlowType zeroFlow=FlowType::Zero();
        //        typedef QuadratureDynamic<1,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadratureDynamicType;
        //        typedef QuadPowDynamic<SplineBase<dim,corder>::pOrder,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadPowDynamicType;
        
        enum MeshLocation{outsideMesh=-1, insideMesh=0, onMeshBoundary=1, onRegionBoundary=2};
        
        
    };
    
    template <int dim>
    struct TypeTraits<DislocationNetwork<dim> > :
    public DislocationNetworkTraitsBase <dim>{};
    
    template <int dim>
    struct TypeTraits<DislocationLoop<dim> > :
    public DislocationNetworkTraitsBase <dim>{};
    
    template <int dim>
    struct TypeTraits<DislocationLoopNode<dim> > :
    public DislocationNetworkTraitsBase <dim>{};
    
    template <int dim>
    struct TypeTraits<DislocationLoopLink<dim> > :
    public DislocationNetworkTraitsBase <dim>{};
    
    template <int dim>
    struct TypeTraits<DislocationNode<dim> > :
    public DislocationNetworkTraitsBase <dim>{};
    
    template <int dim>
    struct TypeTraits<DislocationSegment<dim> > :
    public DislocationNetworkTraitsBase <dim>{};
    
} // namespace model
#endif
