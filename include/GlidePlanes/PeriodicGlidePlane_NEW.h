/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_H_
#define model_PeriodicGlidePlane_H_


#include <memory>
#include <string>
#include <list>

#include <GlidePlane.h>
#include <GlidePlaneFactory.h>
//#include <PeriodicGlidePlaneFactory.h>
#include <TerminalColors.h>
#include <DiophantineSolver.h>
#include <SharedPtrFactories.h>
#include <GramSchmidt.h>
#include <Polygon2D.h>


namespace model
{

template<int dim>
struct PeriodicGlidePlaneFactory;

template<int dim>
struct PeriodicPlanePatch;

template<int dim>
class PeriodicPlaneNode;

template<int dim>
struct PeriodicPatchBoundary;

template<int dim>
class PeriodicGlidePlane;

template<int dim>
struct PeriodicPlaneEdge
{
    
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    const PeriodicPlanePatch<dim>* const patch;
    const std::shared_ptr<PeriodicPlaneNode<dim>> source;
    const std::shared_ptr<PeriodicPlaneNode<dim>>   sink;
    const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection;
    const short int edgeID;
    const VectorDim deltaShift;
    const bool isPeriodicEdge;
    const VectorLowerDim outNormal;
    
    PeriodicPlaneEdge<dim>* next;
    PeriodicPlaneEdge<dim>* prev;
    PeriodicPlaneEdge<dim>* twin;
        
    PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                      const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                      const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                      const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in,
                      const short int&)
    ;
    ~PeriodicPlaneEdge();
    
    static bool periodicFacesInSameGrain(const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in);
    
    std::string tag() const;
    PeriodicPlaneEdge<dim>* nextUntwinned() const;
    PeriodicPlaneEdge<dim>* prevUntwinned() const;
    
    template<typename T>
    std::vector<std::pair<VectorLowerDim,const T* const>> trimPolygon(const std::vector<std::pair<VectorLowerDim,const T* const>>& polyPoints) const
    {
//        // Find if source and sink are inside the original polygon
//        std::vector<VectorLowerDim> windVector;
//        for(const auto& pt : polyPoints)
//        {
//            windVector.push_back(pt.first);
//        }
//        const bool sourceInside(Polygon2D::windingNumber(*source,windVector)!=0);
//        const bool   sinkInside(Polygon2D::windingNumber(  *sink,windVector)!=0);
//
////        NEED TO KEPP EDGE END POINT OF CONTAINED IN POLYGON?
//
////        plane.signedDistanceTo may be better than isAbove/isBelow
        
        std::vector<std::pair<VectorLowerDim,const T* const>> trimmedPolygon;
        const Plane<dim-1> plane(0.5*(*source+*sink),outNormal);
        for (size_t k = 0; k < polyPoints.size(); ++k)
        {
            const auto startPoint(polyPoints[k]);
            const auto endPoint(k == polyPoints.size() - 1 ? polyPoints[0] : polyPoints[k + 1]);
            
            const bool startPointBelowPlane(plane.isAbove(startPoint.first));
            const bool endPointBelowPlane(plane.isAbove(endPoint.first));
            
            if(startPointBelowPlane && endPointBelowPlane)
            {// no crossing, add endPoint
                trimmedPolygon.push_back(endPoint);
            }
            else if(startPointBelowPlane && !endPointBelowPlane)
            {//
                FiniteLineSegment<dim-1> fls(startPoint.first,endPoint.first);
                PlaneSegmentIntersection<dim-1> psi(plane,fls);
                if(psi.type==PlaneSegmentIntersection<dim-1>::INCIDENT)
                {
                    trimmedPolygon.emplace_back(0.5*(psi.x0+psi.x1),nullptr);
                }
                else if(psi.type==PlaneSegmentIntersection<dim-1>::COINCIDENT)
                {// end point may be on the edge, keep it
                    trimmedPolygon.push_back(endPoint);
                }
                else
                {
                    throw std::runtime_error("trimPolygon: intersection should be INCIDENT or COINCIDENT");
                }
            }
            else if(!startPointBelowPlane && endPointBelowPlane)
            {//
                FiniteLineSegment<dim-1> fls(startPoint.first,endPoint.first);
                PlaneSegmentIntersection<dim-1> psi(plane,fls);
                if(psi.type==PlaneSegmentIntersection<dim-1>::INCIDENT)
                {
                    trimmedPolygon.emplace_back(0.5*(psi.x0+psi.x1),nullptr);
                    trimmedPolygon.push_back(endPoint);
                }
                else if(psi.type==PlaneSegmentIntersection<dim-1>::COINCIDENT)
                {// end point may be on the edge, keep it
                    trimmedPolygon.push_back(endPoint);
                }
                else
                {
                    throw std::runtime_error("trimPolygon: intersection should be INCIDENT or COINCIDENT");
                }
            }
            else
            {// both points are above the plane, discard them
                
            }
        }
        return trimmedPolygon;
    }

};

template<int dim>
struct NodalConnectivity
{
    typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
    
    PeriodicPlaneEdgeType* inEdge;
    PeriodicPlaneEdgeType* outEdge;
    
    NodalConnectivity() :
    /* init */ inEdge(nullptr)
    /* init */,outEdge(nullptr)
    {
        
    }
};

template<int dim>
class PeriodicPlaneNode : public StaticID<PeriodicPlaneNode<dim>>
/*                     */,public Eigen::Matrix<double,dim-1,1>
{
    
public:
    
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    typedef PeriodicPlaneNode<dim> PeriodicPlaneNodeType;
    typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
    typedef NodalConnectivity<dim> NodalConnectivityType;
    typedef std::map<size_t,NodalConnectivity<dim>> NodalConnectivityContainerType;
    typedef std::set<PeriodicPlaneEdgeType*> InOutEdgeContainerType;
    typedef VectorLowerDim KeyType;
    typedef CompareVectorsByComponent<double,dim-1,float> CompareType;
    
private:
    
    NodalConnectivityContainerType _patchConnectivities;
    NodalConnectivityContainerType _neighborConnectivities;
    
public:
    
    PeriodicPlaneNode(PeriodicPatchBoundary<dim>* const,const VectorLowerDim& pos);
    void addLink(PeriodicPlaneEdgeType* const link);
    void removeLink(PeriodicPlaneEdgeType* const link);
    const NodalConnectivityContainerType& patchConnectivities() const;
    NodalConnectivityContainerType& patchConnectivities();
    const NodalConnectivityContainerType& neighborConnectivities() const;
    NodalConnectivityContainerType& neighborConnectivities();
    InOutEdgeContainerType inEdges() const;
    InOutEdgeContainerType outEdges() const;
};

template<int dim>
struct PeriodicPlanePatch : public StaticID<PeriodicPlanePatch<dim>>
/*                      */, private std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>>
{
    
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef VectorDim KeyType;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    typedef CompareVectorsByComponent<double,dim,float> CompareType;
    typedef std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>> PeriodicPlaneEdgeContainerType;
    
    PeriodicPatchBoundary<dim>* const patchBoundary;
    const VectorDim shift;
    const std::shared_ptr<GlidePlane<dim>> glidePlane;
    
    PeriodicPlanePatch(PeriodicPatchBoundary<dim>* const,const VectorDim&);
    ~PeriodicPlanePatch();
    static bool isRightHandedBoundary(const BoundingMeshSegments<dim>&, const Plane<dim>&);
    void addMeshIntersections(const BoundingMeshSegments<dim>&);
    const PeriodicPlaneEdgeContainerType& edges() const;
    std::shared_ptr<PeriodicPlaneEdge<dim>> edgeOnFaces(const std::set<const PlanarMeshFace<dim>*>& faces) const;
    int contains(const VectorLowerDim& test) const;
};

template<int dim>
struct PeriodicPatchBoundary : public StaticID<PeriodicPatchBoundary<dim>>
/*                          */,public KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlaneNode<dim>,typename PeriodicPlaneNode<dim>::CompareType>
/*                          */,private std::set<const PeriodicPlaneEdge<dim>*>
/*                          */,public KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlanePatch<dim>,typename PeriodicPlanePatch<dim>::CompareType>

{
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    typedef KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlaneNode<dim>,typename PeriodicPlaneNode<dim>::CompareType> NodeCointainerType;
    typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
    typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
    typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
    typedef KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlanePatch<dim>,typename PeriodicPlanePatch<dim>::CompareType> PatchContainerType;
    
    BoundariesContainerType _outerBoundaries;
    BoundariesContainerType _innerBoundaries;
    GlidePlaneFactory<dim>& glidePlaneFactory;
    const std::shared_ptr<GlidePlane<dim>> referencePlane;

    PeriodicPatchBoundary(GlidePlaneFactory<dim>&, const std::shared_ptr<GlidePlane<dim>>& referencePlane_in);
    void createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy);
    void updateBoundaries();
    BoundariesContainerType& outerBoundaries();
    const BoundariesContainerType& outerBoundaries() const;
    BoundariesContainerType& innerBoundaries();
    const BoundariesContainerType& innerBoundaries() const;
    const UntwinnedEdgeContainerType& untwinnedEdges() const;
    UntwinnedEdgeContainerType& untwinnedEdges();
    void addUntwinnedEdge(const PeriodicPlaneEdge<dim>* link);
    void removeUntwinnedEdge(const PeriodicPlaneEdge<dim>* link);
    NodeCointainerType& nodes();
    const NodeCointainerType& nodes() const;
    int isInsideOuterBoundary(const VectorLowerDim& test);
    bool isCompact() const;
    static VectorDim rightHandedNormal(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
    double rightHandedArea(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
    bool isRightHandedBoundary(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
    VectorLowerDim getLocalPosition(const VectorDim& point,const VectorDim& shift) const;
    VectorDim getGlobalPosition(const VectorLowerDim& point) const;
    std::shared_ptr<PeriodicPlaneNode<dim> > getSharedNode(const VectorDim& pointDim,const VectorDim& shift);
    std::shared_ptr<PeriodicPlanePatch<dim>> getPatch(const VectorDim& shift);
    GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift);
    std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift);
    const PatchContainerType& patches() const;
    PatchContainerType& patches();
};

template<int dim>
class PeriodicGlidePlane : public PeriodicPatchBoundary<dim>
{
        
public:
    
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<long int, dim, 1> VectorDimI;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    //        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> PatchContainerType;
    typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
    typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
    typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
    typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
    typedef GlidePlane<dim> GlidePlaneType;
    typedef typename GlidePlaneType::KeyType KeyType;
    
    GlidePlaneFactory<dim>& glidePlaneFactory;
    PeriodicGlidePlaneFactoryType& periodicGlidePlaneFactory;
    const std::shared_ptr<GlidePlane<dim>> referencePlane;
    
    
    PeriodicGlidePlane(PeriodicGlidePlaneFactoryType* const pgpf,const KeyType& key_in);
    void print();

    
    template<typename T>
    std::vector<std::tuple<VectorLowerDim,VectorDim,std::pair<short int,short  int>,std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t ,const T* const>> polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints,const bool& intersectInternalNodes);
    template<typename T>
    std::vector<std::tuple<VectorLowerDim,VectorDim,std::pair<short int,short  int>,std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t ,const T* const>> polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim,const T* const>>& polyPoints,const bool& intersectInternalNodes);
    template<typename T>
    std::vector<std::tuple<VectorLowerDim,VectorDim,std::pair<short int,short  int>,std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t ,const T* const>> polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints,const bool& intersectInternalNodes,std::set<const PlanarMeshFace<dim>*>& usedFaces);
    template<typename T>
    std::vector<std::tuple<VectorLowerDim,VectorDim,std::pair<short int,short  int>,std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t ,const T* const>> polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim,const T* const>>& polyPoints,const bool& intersectInternalNodes,std::set<const PlanarMeshFace<dim>*>& usedFaces);
    
    std::pair<VectorDim,bool> findPatch(const VectorLowerDim&,const VectorDim&);
    std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> filledPatches(const std::vector<VectorDim>& patchShifts);
};

template<int dim>
struct PeriodicPlanePatchIO
{
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    size_t glidePlaneID;
    size_t patchBoundaryID;
    VectorDim shift;
    
    PeriodicPlanePatchIO();
    PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch);
    PeriodicPlanePatchIO(std::stringstream& ss);
};

template <int dim,class T>
T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds)
{
    os  << ds.glidePlaneID<<"\t"
    /**/<< ds.patchBoundaryID<<"\t"
    /**/<< std::setprecision(15)<<std::scientific<<ds.shift.transpose();
    return os;
}

template<int dim>
template<typename T>
std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints, const bool& intersectInternalNodes)
{
    std::set<const PlanarMeshFace<dim>*> usedFaces;
    return polygonPatchIntersection(polyPoints,intersectInternalNodes,usedFaces);
}

template<int dim>
template<typename T>
std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints, const bool& intersectInternalNodes,std::set<const PlanarMeshFace<dim>*>& usedFaces)
{
    std::vector<std::pair<VectorLowerDim,const T* const>> lowerPolyPoints;
    for(const auto& point : polyPoints)
    {
        assert(this->referencePlane->contains(point.first) && "reference plane does not cointain point");
        lowerPolyPoints.push_back(std::make_pair(referencePlane->localPosition(point.first-VectorDim::Zero()),point.second));
    }
    return polygonPatchIntersection(lowerPolyPoints,intersectInternalNodes,usedFaces);
}

template<int dim>
template<typename T>
std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim, const T* const>>& polyPoints, const bool& intersectInternalNodes)
{
    std::set<const PlanarMeshFace<dim>*> usedFaces;
    return polygonPatchIntersection(polyPoints,intersectInternalNodes,usedFaces);
}

template<int dim>
template<typename T>
std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim, const T* const>>& polyPoints, const bool& intersectInternalNodes,std::set<const PlanarMeshFace<dim>*>& usedFaces)
{
    std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>, std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>> ,size_t,const T* const>> ppiPPedges; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
    // 2dPos,shift,std::pair<edgeIDs>,edgeMap,edgesStillRemaining,LoopNodeType*
    //This is done to take care of all the edges that the two nodes are intersecting, size_t indicates the remaining intersection
    if (polyPoints.size())
    {
        PeriodicPatchBoundary<dim> patchBoundary(glidePlaneFactory, referencePlane);
//        auto lastPatch(patchBoundary.getPatch(findPatch(polyPoints[0].first, VectorDim::Zero()).first));
        std::shared_ptr<PeriodicPlanePatch<dim>> lastPatch(nullptr);
        int startNodeID(-1);
        for (size_t k = 0; k < polyPoints.size(); ++k)
        {
            auto patchFound(findPatch(polyPoints[k].first, VectorDim::Zero()));
            if(patchFound.second)
            {
                lastPatch=patchBoundary.getPatch(patchFound.first);
                startNodeID=k;
                break;
            }
        }
        
        if(lastPatch)
        {
            
//            for(const auto& trimEdge : lastPatch->edges())
//            {
//                if(!trimEdge->isPeriodicEdge)
//                {
//                    polyPoints=trimEdge->trimPolygon(polyPoints);
//                }
//            }
            
            
            // Find other points
            for (size_t k = 0; k < polyPoints.size(); ++k)
            {
                std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short int>,const T* const>> ppi; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
                std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>> edgeMaps; //This will populate the set of edges that this link has crossed along with the patchShift
                std::vector<std::pair<short int,short  int>> edgeVector; //This will just give the number of edges intersected
                std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> currentPatches;
                //            std::cout<<"k point ="<<k<<std::endl;
                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                const auto startPoint(polyPoints[startNodeID+k]);
                const auto endPoint(k == polyPoints.size() - 1 ? polyPoints[startNodeID] : polyPoints[startNodeID+k + 1]);
                if(startPoint.second->periodicPlanePatch() != endPoint.second->periodicPlanePatch() || intersectInternalNodes)
                {
                    while (true)
                    {
    //                    std::map<const PeriodicPlanePatch<dim> *, std::vector<std::pair<const VectorLowerDim, const PeriodicPlaneEdge<dim> *const>>> crossdEdges;
//                        std::vector<std::pair<const VectorLowerDim, const PeriodicPlaneEdge<dim> *const>> crossdEdges; // TODO this should be a container sorted by the scalae intersection paramenter along startPoint->endPoint (ssd.t)

                        std::map<float,std::pair<VectorLowerDim,std::set<const PeriodicPlaneEdge<dim>*>>> crossdEdges;
                        
                        for (const auto &bndEdge : patchBoundary.untwinnedEdges())
                        { // loop over outer boundaries and holes
                            
                            SegmentSegmentDistance<dim - 1> ssd(startPoint.first, endPoint.first, *bndEdge->source, *bndEdge->sink);
                            if (ssd.isIntersecting(FLT_EPSILON))
                            { // intersection with current boundary found
                                PlaneSegmentIntersection<dim-1> psi(0.5*(*bndEdge->source+*bndEdge->sink),bndEdge->outNormal,startPoint.first,endPoint.first);
                                if(   psi.type==PlaneSegmentIntersection<dim-1>::INCIDENT
                                   && !bndEdge->isPeriodicEdge)
                                {// try to trim polygon
                                    for(const auto& face : bndEdge->meshIntersection->faces)
                                    {
                                        const bool isTrimmingFace(face->periodicFacePair.second==nullptr? true : face->regionIDs!=face->periodicFacePair.second->regionIDs);
                                        if(isTrimmingFace)
                                        {// not a periodic face
                                            if(usedFaces.find(face)==usedFaces.end())
                                            {// face not used yet
                                                usedFaces.insert(face);
                                                return polygonPatchIntersection(bndEdge->trimPolygon(polyPoints),intersectInternalNodes,usedFaces);
                                            }
                                        }
                                    }
                                }
                                
                                // Store intersection(s) ordered along loop segment
                                const auto interSegment(ssd.intersectionSegment());
                                for(const auto& tup : interSegment)
                                {
                                    const auto crossIter(crossdEdges.find(std::get<1>(tup)));
                                    if(crossIter==crossdEdges.end())
                                    {
                                        std::set<const PeriodicPlaneEdge<dim>*> tempSet;
                                        tempSet.insert(bndEdge);
                                        crossdEdges.emplace(std::get<1>(tup),std::get<0>(tup),tempSet);
                                    }
                                    else
                                    {
                                        crossIter->second.second.insert(bndEdge);
                                    }
                                }
                            }
                        }
                        
                        if (crossdEdges.size() == 0)
                        { // No untwinned edges have been crossed, break and move to next pair of polyPoints
                            ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                            break; // while loop
                            
                            THERE IS A PROBLEM HERE BACAUSE THIS IS THE ONLY BREAK FROM THE MAIN WHILE LOOP
                            HOWEVER, FOR A TRMMING EDDE WE KEEP FINDING INTERSECTIONS AND WE NEVER GET HERE
                            
                        }
                        else
                        {
                            // START OF NEW IMPLEMENTATION
                            switch (crossdEdges.size())
                            {
                
                                    CONTINUE HERE
                                    
                                    for(const auto& crossdEdgePair : crossdEdges)
                                    {
//                                        if(crossdEdgePair.second.second.size())
//                                        {
//                                            const auto& edgePatch(crossdEdgePair.second.second.begin()->patch);
//                                            for(const auto& otherEdge : crossdEdgePair.second.second())
//                                            {
//                                                if(otherEdge->paatch!=edgePatch)
//                                                {
//                                                    throw std::runtime_error("Two different patches are found for the same edge intersection");
//                                                }
//                                            }
//                                        }
                                    }
                                    
                                case 1:
                                {// only one intersection point (possibly one or two edges)
                                    
                                    
                                    const auto& edge(*crossdEdges.begin());
                                    const auto& edgePatchShift(edge.second->patch->shift);

                                    ppi.emplace_back(edge.first, edgePatchShift, std::make_pair(edge.second->edgeID,-1), nullptr);
                                    edgeVector.push_back(std::make_pair(edge.second->edgeID,-1));
                                    
                                    const auto edgeSetIter1(edgeMaps.find(edgePatchShift));
                                    if (edgeSetIter1 == edgeMaps.end())
                                    {
                                        std::set<std::pair<short int,short  int>> edgeSet;
                                        edgeSet.insert(std::make_pair(edge.second->edgeID,-1));
                                        edgeMaps.emplace(edgePatchShift, edgeSet);
                                    }
                                    else
                                    {
                                        edgeSetIter1->second.insert(std::make_pair(edge.second->edgeID,-1));
                                    }
                                    
                                    if(edge.second->isPeriodicEdge)
                                    {
                                        const VectorDim localShift(edgePatchShift + edge.second->deltaShift);
                                        // Insert patches corresponding to the individual intersections as well
                                        lastPatch = patchBoundary.patches().getFromKey(localShift);
                                        currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                        assert(edge.second->twin);
                                        ppi.emplace_back(edge.first, edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1), nullptr);
                                        edgeVector.push_back(std::make_pair(edge.second->twin->edgeID,-1));
                                        
                                        const auto edgeSetIter2(edgeMaps.find(edge.second->twin->patch->shift));
                                        if (edgeSetIter2 == edgeMaps.end())
                                        {
                                            std::set<std::pair<short int,short  int>> edgeSet;
                                            edgeSet.insert(std::make_pair(edge.second->twin->edgeID,-1));
                                            edgeMaps.emplace(edge.second->twin->patch->shift, edgeSet);
                                        }
                                        else
                                        {
                                            edgeSetIter2->second.insert(std::make_pair(edge.second->twin->edgeID,-1));
                                        }
                                    }
                                    else
                                    {
    //                                    Line<dim-1> line(edge.first,(*edge.second->sink)-(*edge.second->source));
    //                                    std::pair<double,double> lineBounds(std::make_pair(0.0,0.0));
    //
    //                                    for (size_t j = k+1; j < polyPoints.size(); ++k)
    //                                    {// Pick subsequent loop segments until you find another intersection with the line
    //                                        const auto startPointTemp(polyPoints[startNodeID+j]);
    //                                        const auto endPointTemp(j == polyPoints.size() - 1 ? polyPoints[startNodeID] : polyPoints[startNodeID+j + 1]);
    //
    //                                        const auto chord(endPointTemp.first-startPointTemp.first);
    //                                        const double chordNorm(chord.norm());
    //                                        if(chordNorm>FLT_EPSILON)
    //                                        {
    //                                            LineLineIntersection<dim-1> lli(line.P,line.d,startPointTemp.first,chord/chordNorm);
    //                                            if(   lli.type==LineLineIntersection<dim-1>::INCIDENT
    //                                               && lli.u1/chordNorm > 0.0-FLT_EPSILON // intersection inside finite line segment
    //                                               && lli.u1/chordNorm < 1.0+FLT_EPSILON)// intersection inside finite line segment
    //                                            {
    //                                                lineBounds=std::make_pair(std::min(0.0,lli.u0),std::max(0.0,lli.u0));
    //                                                break;
    //                                            }
    //                                        }
    //                                    }
    //                                    if(lineBounds.first==lineBounds.second)
    //                                    {
    //                                        throw std::runtime_error("Cound not determine lineBounds.");
    //                                    }
    //
    //                                    // Check if polyPoints polygon
    //
    //                                    std::set<const PeriodicPlaneEdge<dim>*> possibleEdges;
    //                                    for (const auto &bndEdge : patchBoundary.untwinnedEdges()) // instead of this for loop we could use only edge.second->nextUntwinned() and edge.second->prevUntwinned()
    //                                    { // loop over outer boundaries and holes
    //                                        const auto chord(*bndEdge->sink-*bndEdge->source);
    //                                        const double chordNorm(chord.norm());
    //                                        if(chordNorm>FLT_EPSILON)
    //                                        {
    //                                            LineLineIntersection<dim-1> lli(line.P,line.d,*bndEdge->source,chord/chordNorm);
    //                                            if(   lli.type==LineLineIntersection<dim-1>::INCIDENT
    //                                               && lli.u0 > lineBounds.first // intersection inside lineBounds
    //                                               && lli.u0 < lineBounds.second // intersection inside lineBounds
    //                                               && lli.u1/chordNorm > 0.0-FLT_EPSILON // intersection inside finite line segment
    //                                               && lli.u1/chordNorm < 1.0+FLT_EPSILON)// intersection inside finite line segment
    //                                            {
    //                                                possibleEdges.insert(bndEdge);
    //                                            }
    //                                        }
    //                                    }
    //
    //                                    if(possibleEdges.size()==1)
    //                                    {
    //
    //                                    }
    //                                    else
    //                                    {
    //                                        std::cout<<"possibleEdges.size()="<<possibleEdges.size()<<std::endl;
    //                                        throw std::runtime_error("Expected exactly one possible edge");
    //                                    }
                                        
                                    }
                                    
                                    
                                    break;
                                }
                                    
                                case 2:
                                {
                                    const auto& edge0(*crossdEdges.begin());
                                    const auto& edge1(*crossdEdges.rbegin());
                                    if(   (edge0.first-edge1.first).norm()>FLT_EPSILON // 2d intersection not coincident
                                       || edge0.second==edge1.second // same edge twice
                                       || edge0.second->patch!=edge1.second->patch // edges in different patches
                                       )
                                    {
                                        throw std::runtime_error("polygonPatchIntersection: two intersections are inconsistent.");
                                    }
                                    const auto& edgePatchShift(edge0.second->patch->shift);
                                    
                                    ppi.emplace_back(edge0.first, edgePatchShift, std::make_pair(edge0.second->edgeID, edge1.second->edgeID), nullptr);
                                    edgeVector.push_back(std::make_pair(edge0.second->edgeID, edge1.second->edgeID));

                                    const auto edgeSetIter1 (edgeMaps.find(edgePatchShift));
                                    if (edgeSetIter1==edgeMaps.end())
                                    {
                                        std::set<std::pair<short int,short  int>> edgeSet;
                                        edgeSet.insert(std::make_pair(edge0.second->edgeID, edge1.second->edgeID));
                                        edgeMaps.emplace(edgePatchShift,edgeSet);
                                    }
                                    else
                                    {
                                        edgeSetIter1->second.insert(std::make_pair(edge0.second->edgeID, edge1.second->edgeID));
                                    }
                                    
                                    //here three patches are needed to be inserted... two at periodically opposite faces and one at diagonally opposite faces
                                    lastPatch = patchBoundary.patches().getFromKey(edgePatchShift + edge0.second->deltaShift);
                                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                    lastPatch = patchBoundary.patches().getFromKey(edgePatchShift + edge1.second->deltaShift);
                                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                    const VectorDim localShift(edgePatchShift + edge0.second->deltaShift+edge1.second->deltaShift);
                                    lastPatch = patchBoundary.patches().getFromKey(localShift);
                                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                    assert(edge0.second->twin); //Now the twin for both the edges must exist
                                    assert(edge1.second->twin);
                                    assert((edge0.second->twin!=edge1.second->twin) && "Two intersection must be different after twinned link");
                                    //Insert a new node at the twinned edge of the current pair of edges
                                    //Need to determine the correct edgeIDs
                                    short int beginTwinEdgeID(-1);
                                    short rbeginTwinEdgeID(-1);
                                    //We need to get the edgeID corresponding to the diagonally opposite patch
                                    if (edge0.second==edge1.second->prev)
                                    {
                                        assert(edge0.second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");
                                        assert(edge1.second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
                                        beginTwinEdgeID=edge0.second->twin->prev->twin->edgeID;
                                        rbeginTwinEdgeID=edge1.second->twin->next->twin->edgeID;
                                    }
                                    else if (edge0.second==edge1.second->next)
                                    {
                                        assert(edge0.second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
                                        assert(edge1.second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");
                                        
                                        beginTwinEdgeID=edge0.second->twin->next->twin->edgeID;
                                        rbeginTwinEdgeID=edge1.second->twin->prev->twin->edgeID;
                                    }
                                    else
                                    {
                                        assert(false && "Edges must be consecutive");
                                    }
                                    
                                    assert(beginTwinEdgeID!=-1 && "A twin ID must exist");
                                    assert(rbeginTwinEdgeID!=-1 && "A twin ID must exist");
                                    
                                    ppi.emplace_back(edge0.first, localShift, std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID), nullptr);
                                    edgeVector.push_back(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                                    
                                    const auto edgeSetIter2 (edgeMaps.find(localShift));
                                    if (edgeSetIter2==edgeMaps.end())
                                    {
                                        std::set<std::pair<short int,short  int>> edgeSet;
                                        edgeSet.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                                        edgeMaps.emplace(localShift,edgeSet);
                                    }
                                    else
                                    {
                                        edgeSetIter2->second.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                                    }
                                    
                                    break;
                                }
                                    

                                    
                                default:
                                {
                                    std::cout<<"crossdEdges.size()="<<crossdEdges.size()<<std::endl;
                                    throw std::runtime_error("polygonPatchIntersection: at most two intersections are possible.");
                                    break;
                                }
                            }
                            
                            // OLD IMPLEMENTATION
    //                        for (const auto &pair : crossdEdges) //HOW IS THIS IN THE RIGHT ORDER???
    //                        { // loop over crossed patches
    //                            assert((pair.second.size()>0 && pair.second.size()<=2) && "At max two periodic plane edges can be intersected (2D geometry constraint)");
    //                            if (pair.second.size()==2 && (pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)
    //                            {// Two intersectoins at the same position (intersectoin at corner)
    //                                // pair.second.size() cannot be 4 even for diagonally opposite ends since we are looping over the untwineed edges for the intersection.
    //                                //Need to check that two subsequen links are at the same 2D position
    //
    //                                //The two intersection points are the same... Insert at diagonally opposite end
    //
    //                                assert((pair.second.begin()->second!=pair.second.rbegin()->second) && "Two intersection must be different");
    //                                ppi.emplace_back(pair.second.begin()->first, pair.first->shift, std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID), nullptr);
    //                                const auto edgeSetIter1 (edgeMaps.find(pair.first->shift));
    //                                edgeVector.push_back(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
    //
    //                                if (edgeSetIter1==edgeMaps.end())
    //                                {
    //                                    std::set<std::pair<short int,short  int>> edgeSet;
    //                                    edgeSet.insert(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
    //                                    edgeMaps.emplace(pair.first->shift,edgeSet);
    //                                }
    //                                else
    //                                {
    //                                    edgeSetIter1->second.insert(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
    //                                }
    //
    //                                //here three patches are needed to be inserted... two at periodically opposite faces and one at diagonally opposite faces
    //                                lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.begin()->second->deltaShift);
    //                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
    //                                lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.rbegin()->second->deltaShift);
    //                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
    //                                const VectorDim localShift(pair.first->shift + pair.second.begin()->second->deltaShift+pair.second.rbegin()->second->deltaShift);
    //                                lastPatch = patchBoundary.patches().getFromKey(localShift);
    //                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
    //                                assert(pair.second.begin()->second->twin); //Now the twin for both the edges must exist
    //                                assert(pair.second.rbegin()->second->twin);
    //                                assert((pair.second.begin()->second->twin!=pair.second.rbegin()->second->twin) && "Two intersection must be different after twinned link");
    //                                //Insert a new node at the twinned edge of the current pair of edges
    //                                //Need to determine the correct edgeIDs
    //                                short int beginTwinEdgeID(-1);
    //                                short rbeginTwinEdgeID(-1);
    //                                //We need to get the edgeID corresponding to the diagonally opposite patch
    //                                if (pair.second.begin()->second==pair.second.rbegin()->second->prev)
    //                                {
    //                                    assert(pair.second.begin()->second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");
    //                                    assert(pair.second.rbegin()->second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
    //                                    beginTwinEdgeID=pair.second.begin()->second->twin->prev->twin->edgeID;
    //                                    rbeginTwinEdgeID=pair.second.rbegin()->second->twin->next->twin->edgeID;
    //                                }
    //                                else if (pair.second.begin()->second==pair.second.rbegin()->second->next)
    //                                {
    //                                    assert(pair.second.begin()->second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
    //                                    assert(pair.second.rbegin()->second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");
    //
    //                                    beginTwinEdgeID=pair.second.begin()->second->twin->next->twin->edgeID;
    //                                    rbeginTwinEdgeID=pair.second.rbegin()->second->twin->prev->twin->edgeID;
    //                                }
    //                                else
    //                                {
    //                                    assert(false && "Edges must be consecutive");
    //                                }
    //
    //                                assert(beginTwinEdgeID!=-1 && "A twin ID must exist");
    //                                assert(rbeginTwinEdgeID!=-1 && "A twin ID must exist");
    //
    //                                ppi.emplace_back(pair.second.begin()->first, localShift, std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID), nullptr);
    //                                edgeVector.push_back(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
    //
    //                                const auto edgeSetIter2 (edgeMaps.find(localShift));
    //                                if (edgeSetIter2==edgeMaps.end())
    //                                {
    //                                    std::set<std::pair<short int,short  int>> edgeSet;
    //                                    edgeSet.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
    //                                    edgeMaps.emplace(localShift,edgeSet);
    //                                }
    //                                else
    //                                {
    //                                    edgeSetIter2->second.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
    //                                }
    //                            }
    //                            else
    //                            {// Intersections not a corner
    //                                for (const auto &edge : pair.second) HOW IS THIS IN THE RIGHT ORDER???
    //                                { // loop over crossed edges. edge.first=2D position, edge.second= edge ptr
    //                                    ppi.emplace_back(edge.first, pair.first->shift, std::make_pair(edge.second->edgeID,-1), nullptr);
    //                                    edgeVector.push_back(std::make_pair(edge.second->edgeID,-1));
    //
    //                                    const auto edgeSetIter1(edgeMaps.find(pair.first->shift));
    //                                    if (edgeSetIter1 == edgeMaps.end())
    //                                    {
    //                                        std::set<std::pair<short int,short  int>> edgeSet;
    //                                        edgeSet.insert(std::make_pair(edge.second->edgeID,-1));
    //                                        edgeMaps.emplace(pair.first->shift, edgeSet);
    //                                    }
    //                                    else
    //                                    {
    //                                        edgeSetIter1->second.insert(std::make_pair(edge.second->edgeID,-1));
    //                                    }
    //                                    const VectorDim localShift(pair.first->shift + edge.second->deltaShift);
    //                                    // Insert patches corresponding to the individual intersections as well
    //                                    lastPatch = patchBoundary.patches().getFromKey(localShift);
    //                                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
    //                                    assert(edge.second->twin);
    //                                    ppi.emplace_back(edge.first, edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1), nullptr);
    //                                    edgeVector.push_back(std::make_pair(edge.second->twin->edgeID,-1));
    //
    //                                    const auto edgeSetIter2(edgeMaps.find(edge.second->twin->patch->shift));
    //                                    if (edgeSetIter2 == edgeMaps.end())
    //                                    {
    //                                        std::set<std::pair<short int,short  int>> edgeSet;
    //                                        edgeSet.insert(std::make_pair(edge.second->twin->edgeID,-1));
    //                                        edgeMaps.emplace(edge.second->twin->patch->shift, edgeSet);
    //                                    }
    //                                    else
    //                                    {
    //                                        edgeSetIter2->second.insert(std::make_pair(edge.second->twin->edgeID,-1));
    //                                    }
    //                                }
    //                            }
    //                        }
                            if (lastPatch->contains(endPoint.first))
                            {
                                ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                                 break;
                            }
                        }
                    }
                }
                else
                {
                    ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                }
                
                //Here populate the ppiPPedges
                size_t crossCount=0;
                for (const auto& ppiIter : ppi )
                {
                    size_t edgesStillRemaining(0);
                    if (std::get<3>(ppiIter) == nullptr)
                    {
                        //That mean intersection with the edge is here
                        crossCount++;
                        edgesStillRemaining = edgeVector.size() - crossCount;
                    }
                    ppiPPedges.emplace_back(std::get<0>(ppiIter), std::get<1>(ppiIter), std::get<2>(ppiIter), edgeMaps, edgesStillRemaining, std::get<3>(ppiIter));
                }
                if (edgeVector.size() != 0)
                {
                    if ((edgeVector.size() - crossCount ) != 0)
                    {
                        std::cout<<"Edge VectorSize "<<edgeVector.size()<<std::endl;
                        std::cout<<"Cross count "<<crossCount<<std::endl;
                    }
                    assert((edgeVector.size() - crossCount) == 0 && "No edges should be remaining after completion of polygonPatchIntersection");
                }
            }
        }

        
        
        
//        auto lastPatch(patchBoundary.getPatch(findPatch(polyPoints[startNodeID].first, VectorDim::Zero()).first));

        // const size_t loopID((*polyPoints.begin()).second->loop()->sID);
        // const size_t runID((*polyPoints.begin()).second->network().simulationParameters.runID);
        
        //                std::ofstream polyFile("poly.txt");
        //                for(const auto& node : polyPoints)
        //                {
        //                    polyFile<<"    "<<node.transpose()<<std::endl;
        //                }
        
        //        int removeMe=0;
        
        //        ppi.emplace_back(polyPoints[0],lastPatch,nullptr);
        

    }
    
    return ppiPPedges;
}

template<typename T1,typename T2, typename T3,typename T4>
void printIntersections(const T1& ppi,const T2& usedPatches,const T3& lastPatch,int& removeMe, const T4& untwinnedEdges)
{
    std::ofstream pointsFile("points_"+std::to_string(removeMe)+".txt");
    for(const auto& tup : ppi)
    {
        //                        if(!node.second.expired())
        //                        {
        pointsFile<<"    "<<std::get<0>(tup).transpose()<<std::endl;
        //                        }
    }
    
    std::ofstream edgesFile("edges_"+std::to_string(removeMe)+".txt");
    for(const auto& patch : usedPatches)
    {
        for(const auto& edge : patch->edges())
        {
            if(untwinnedEdges.find(edge.get())!=untwinnedEdges.end())
            {
                edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 1"<<std::endl;
            }
            else
            {
                edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 0"<<std::endl;
            }
            
        }
    }
    
    std::ofstream lastPatchFile("lastPatch_"+std::to_string(removeMe)+".txt");
    for(const auto& edge : lastPatch->edges())
    {
        lastPatchFile<<lastPatch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<std::endl;
    }
    removeMe++;
}


}
#endif
