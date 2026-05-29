/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationSegment_cpp_
#define model_DislocationSegment_cpp_

#include <DislocationSegment.h>
#include <PlanesIntersection.h>

namespace model
{
    template <int dim>
    DislocationSegment<dim>::DislocationSegment(LoopNetworkType* const net,
                                                                         const std::shared_ptr<NetworkNodeType>& nI,
                                                                         const std::shared_ptr<NetworkNodeType>& nJ) :
    /* init */ NetworkLink<DislocationSegment<dim>>(net,nI,nJ)
    /* init */,SplineSegmentType(this->source->get_P(),this->sink->get_P())
    // /* init */,ConfinedDislocationObjectType(this->source->get_P(),this->sink->get_P())
    /* init */,Burgers(VectorDim::Zero())
    /* init */,BurgersNorm(Burgers.norm())
    /* init */,straight(this->network().ddBase.poly,this->source->get_P(),this->sink->get_P(),Burgers,this->chordLength(),this->unitDirection(),this->network().ddBase.EwaldLength)
    /* init */,_slipSystem(nullptr)
    {
        VerboseDislocationSegment(1,"Constructing DislocationSegment "<<this->tag()<<std::endl);
    }
    
template <int dim>
typename DislocationSegment<dim>::VectorDim DislocationSegment<dim>::climbDirection() const
{
    if(this->isBoundarySegment()
       || this->isGrainBoundarySegment())
    {
        return VectorDim::Zero();
    }
    else
    {
        const VectorDim temp(this->burgers().cross(this->chord()));
        const double tempNorm(temp.norm());
        return tempNorm<FLT_EPSILON? VectorDim::Zero() : (temp/tempNorm).eval();
    }
}

       template <int dim>
    void DislocationSegment<dim>::updateSlipSystem()
    {
//        if (this->grains().size() == 0)
//        {
//             _slipSystem=nullptr;
//        }
        const auto grs(this->grains());
        if (grs.size() == 1)
        {
            std::set<std::shared_ptr<SlipSystem>> ssSet;

            for (const auto &loopLink : this->loopLinks())
            {
                if (loopLink->loop->slipSystem())
                {
                    ssSet.insert(loopLink->loop->slipSystem());
                }
            }
            
            switch (ssSet.size())
            {
                case 0:
                {
                    _slipSystem = nullptr;
                    break;
                }
                    
                case 1:
                {
                    _slipSystem = *ssSet.begin();
                    break;
                }
                    
                default:
                {
                    _slipSystem = nullptr;

                    const auto gps(glidePlanes());
                    if(gps.size()==1)
                    {//Possible planar glissile junction, we'll pick the first slip system that has the same burgers
                        for(const auto& gpss : (*gps.begin())->slipSystems())
                        {
                            if((gpss->s.cartesian()-burgers()).norm()<FLT_EPSILON)
                            {
                                _slipSystem=gpss;
                                break;
                            }
                        }
                    }
                    break;
                }
            }

//            if (ssSet.size() == 1)
//            { // a unique slip system found.
//                _slipSystem = *ssSet.begin();
//            }
//            else
//            { // ssSet.size()==0.
//                _slipSystem = nullptr;
//
//                const auto gps(glidePlanes());
//                if(gps.size()==1)
//                {//Planar glissile junction, we'll pick the first slip system that has the same flow
//                    for(const auto& gpss : (*gps.begin())->slipSystems())
//                    {
//
//                    }
//                }
//                else
//                {
//                    _slipSystem = nullptr;
//                }
//            }
        }
        else
        {
//            assert(false && "FINISH THIS FOR MULTIPLE GRAINS");
            _slipSystem = nullptr;
        }
        
        VerboseDislocationSegment(3, "_slipSystem.s= " << _slipSystem->s.cartesian().transpose() << std::endl;);
        if (_slipSystem)
        {

            VerboseDislocationSegment(3, "_slipSystem.s= " << _slipSystem->s.cartesian().transpose() << std::endl;);
            VerboseDislocationSegment(3, "_slipSystem.n= " << _slipSystem->unitNormal.transpose() << std::endl;);
        }
        else
        {
        }
    }

    /**********************************************************************/
    template <int dim>
    void DislocationSegment<dim>::addLoopLink(LoopLinkType* const pL)
    {
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", adding LoopLink "<<pL->tag()<<std::endl;);
        NetworkLink<DislocationSegment>::addLoopLink(pL); // forward to base class
        if(pL->source->networkNode==this->source)
        {// Update Burgers vector
            Burgers+=pL->flow().cartesian();
        }
        else
        {
            Burgers-=pL->flow().cartesian();
        }
        BurgersNorm=Burgers.norm();
        straight.updateGeometry();
        
        //assert(0 && "Following addGlidePlane is wrong, must consider shifts");
        
        if(this->network().ddBase.isPeriodicDomain)
        {
            const auto periodicPlanePatch(pL->periodicPlanePatch());
            if(periodicPlanePatch)
            {
                // this->confinedObject().addGlidePlane(periodicPlanePatch->glidePlane.get());
               // /* Check for the consistency of the glide planes
                //assert if the confined object glide planes are the same as that of the each node glide planes
                bool glidePlanesConsistency(false);
                for (const auto &Seggps : glidePlanes())
                {
                    glidePlanesConsistency = false;
                    for (const auto &sourceGPs : this->source->glidePlanes())
                    {
                        if (Seggps == sourceGPs)
                        {
                            glidePlanesConsistency = true;
                            break;
                        }
                    }
                }
                // std::cout<<" A " <<glidePlanesConsistency<<std::endl;

                if (glidePlanesConsistency)
                {
                    for (const auto &Seggps : glidePlanes())
                    {
                        glidePlanesConsistency=false;
                        for (const auto &sinkGPs : this->sink->glidePlanes())
                        {
                            if (Seggps==sinkGPs)
                            {
                                glidePlanesConsistency=true;
                                break;
                            }
                        }
                    }
                }
                // std::cout<<" B " <<glidePlanesConsistency<<std::endl;

                assert(glidePlanesConsistency && "Not a consistent definition of glide planes for segments and segment nodes");
              // Check for the consistency of glide planes end here  */
            }
        }
        // else
        // {
        //     this->confinedObject().addGlidePlane(pL->loop->glidePlane.get());
        // }


        updateSlipSystem();
    }
    
    
    
    /**********************************************************************/
    template <int dim>
    void DislocationSegment<dim>::removeLoopLink(LoopLinkType* const pL)
    {
       VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", removing LoopLink "<<pL->tag()<<std::endl;);
        NetworkLink<DislocationSegment>::removeLoopLink(pL);  // forward to base class
        
        if(pL->source->networkNode==this->source)
        {// Modify Burgers vector
            Burgers-=pL->flow().cartesian();
        }
        else
        {
            Burgers+=pL->flow().cartesian();
        }
        BurgersNorm=Burgers.norm();
        straight.updateGeometry(); // update b x t
        
        updateSlipSystem();
        
    }
    
    template <int dim>
    void DislocationSegment<dim>::updateGeometry()
    {
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", updateGeometry "<<std::endl;);
        SplineSegmentType::updateGeometry();
        // this->updateConfinement(this->source->get_P(),this->sink->get_P());
        straight.updateGeometry();
        VerboseDislocationSegment(2,"DislocationSegment "<<this->tag()<<", updateGeometry DONE"<<std::endl;);
    }
    
    template <int dim>
    bool DislocationSegment<dim>::isGlissile() const
    {/*\returns true if ALL the following conditions are met
      * - the segment is confined by only one plane
      * - all loops containing this segment are glissile
      */
        //            bool temp(this->glidePlanes().size()==1 && !hasZeroBurgers() && !isVirtualBoundarySegment());
        bool temp(this->glidePlanes().size()==1);
        if(temp)
        {
            for(const auto& loopLink : this->loopLinks())
            {
                temp=(temp && loopLink->loop->loopType==DislocationLoopIO<dim>::GLISSILELOOP);
            }
        }
        return temp;
    }
    
    template <int dim>
    bool DislocationSegment<dim>::isSessile() const
    {
        return !isGlissile();
    }
    
    /******************************************************************/
    template <int dim>
    void DislocationSegment<dim>::initFromFile(const std::string& fileName)
    {
        //        LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
        //        assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
        //        assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
        quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",false);
        //            assembleWithTangentProjection=TextFileParser(fileName).readScalar<int>("assembleWithTangentProjection",true);
        assert((NetworkLinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
        verboseDislocationSegment=TextFileParser(fileName).readScalar<int>("verboseDislocationSegment",false);
    }
    
    template <int dim>
    typename DislocationSegment<dim>::MeshLocation DislocationSegment<dim>::meshLocation() const
    {/*!\returns the position of *this relative to the bonudary:
      * 1 = inside mesh
      * 2 = on mesh boundary
      */
        if(isBoundarySegment())
        {
            return MeshLocation::onMeshBoundary;
        }
        else if(isGrainBoundarySegment())
        {
            return MeshLocation::onRegionBoundary;
        }
//        else if(isVirtualBoundarySegment())
//        {
//            return MeshLocation::outsideMesh;
//        }
        else
        {
            return MeshLocation::insideMesh;
        }
    }
    
//    template <int dim>
//    bool DislocationSegment<dim>::isVirtualBoundarySegment() const
//    {//!\returns true if all loops of this segment are virtualBoundaryLoops
//        bool temp(true);
//        for(const auto& loopLink : this->loopLinks())
//        {
//            temp= (temp && loopLink->loop->isVirtualBoundaryLoop());
//            if(!temp)
//            {
//                break;
//            }
//        }
//        return temp;
//    }
    
    template <int dim>
    bool DislocationSegment<dim>::hasZeroBurgers() const
    {
        return BurgersNorm<FLT_EPSILON;
    }
    
    template <int dim>
    const std::shared_ptr<SlipSystem>&  DislocationSegment<dim>::slipSystem() const
    {
        return _slipSystem;
    }
    
    template <int dim>
    bool DislocationSegment<dim>::isBoundarySegment() const
    {/*!\returns true if both nodes are boundary nodes, and the midpoint is
      * on the boundary.
      */
        return this->isOnExternalBoundary();
    }
    
    template <int dim>
    bool DislocationSegment<dim>::isGrainBoundarySegment() const
    {
        return this->isOnInternalBoundary();
    }
    
    template <int dim>
    const typename DislocationSegment<dim>::VectorDim& DislocationSegment<dim>::burgers() const
    {
        return Burgers;
    }
    
    template <int dim>
    void DislocationSegment<dim>::addToGlobalAssembly(std::deque<Eigen::Triplet<double> >& kqqT,
                             Eigen::VectorXd& FQ) const
    {/*!\param[in] kqqT the stiffness matrix of the network component
      * \param[in] FQ the force vector of the network component
      */
        if(!hasZeroBurgers())
        {
//            h2posMap.clear();
            std::map<size_t,
            /*    */ std::pair<VectorNcoeff,VectorDim>,
            /*    */ std::less<size_t>
            /*    */ > h2posMap;
            h2posMap.emplace(this->source->networkID(),std::make_pair((VectorNcoeff()<<1.0,0.0).finished(),this->source->get_P()));
            h2posMap.emplace(this->  sink->networkID(),std::make_pair((VectorNcoeff()<<0.0,1.0).finished(),this->  sink->get_P()));

            
//            switch (corder)
//            {
//                case 0:
//                {
//                    h2posMap.emplace(this->source->networkID(),std::make_pair((VectorNcoeff()<<1.0,0.0).finished(),this->source->get_P()));
//                    h2posMap.emplace(this->  sink->networkID(),std::make_pair((VectorNcoeff()<<0.0,1.0).finished(),this->  sink->get_P()));
//                    break;
//                }
//                    
//                default:
//                {
//                    assert(0 && "IMPLEMENT THIS CASE FOR CURVED SEGMENTS");
//                    break;
//                }
//            }
            //        h2posMap=this->hermite2posMap();
            Eigen::Matrix<double, Ndof, Eigen::Dynamic> Mseg(Eigen::Matrix<double, Ndof, Eigen::Dynamic>::Zero(Ncoeff*dim,h2posMap.size()*dim));
//            Mseg.setZero(Ncoeff*dim,h2posMap.size()*dim);
            size_t c=0;
            for(const auto& pair : h2posMap)
            {
                for(int r=0;r<Ncoeff;++r)
                {
                    Mseg.template block<dim,dim>(r*dim,c*dim)=pair.second.first(r)*MatrixDim::Identity();
                }
                c++;
            }
                    

            const Eigen::MatrixXd tempKqq(Mseg.transpose()*this->nodalVelocityMatrix(*this)*Mseg); // Create the temporaty stiffness matrix and push into triplets
            size_t localI=0;
            for(const auto& pairI : h2posMap)
            {
                for(int dI=0;dI<dim;++dI)
                {
                    const size_t globalI=pairI.first*dim+dI;
                    size_t localJ=0;
                    for(const auto& pairJ : h2posMap)
                    {
                        for(int dJ=0;dJ<dim;++dJ)
                        {
                            const size_t globalJ=pairJ.first*dim+dJ;
                            
                            if (std::fabs(tempKqq(localI,localJ))>FLT_EPSILON)
                            {
                                kqqT.emplace_back(globalI,globalJ,tempKqq(localI,localJ));
                            }
                            localJ++;
                        }
                    }
                    localI++;
                }
            }
            
            if(this->quadraturePoints().size())
            {
                const Eigen::VectorXd tempFq(Mseg.transpose()*this->nodalVelocityVector(*this)); // Create temporary force vector and add to global FQ
                localI=0;
                for(const auto& pairI : h2posMap)
                {
                    for(int dI=0;dI<dim;++dI)
                    {
                        const size_t globalI=pairI.first*dim+dI;
                        
                        FQ(globalI)+=tempFq(localI);
                        
                        localI++;
                    }
                }
            }
        }
    }
    
    template <int dim>
    void DislocationSegment<dim>::createQuadraturePoints(const bool& isClimbStep)
    {
        this->create(*this,quadPerLength,isClimbStep);
    }

    template <int dim>
    void DislocationSegment<dim>::updateQuadraturePoints(const bool& isClimbStep)
    {
        this->update(*this,isClimbStep);
    }

    template <int dim>
    int DislocationSegment<dim>::velocityGroup(const double &maxVelocity, const std::set<int> &subcyclingBins) const
    // This represents how many number of steps taken for updating the velocity
    {
        if (this->quadraturePoints().size() == 0)
        {
            return *subcyclingBins.begin();
        }
        else
        {
            if (maxVelocity > FLT_EPSILON)
            {
                const double avgV((0.5 * (this->source->get_V() + this->sink->get_V())).norm());

                if (avgV >= FLT_EPSILON)
                {
                    const double velRat(maxVelocity / avgV);
                    std::map<double, int> temp;
                    for (const auto &sbin : subcyclingBins)
                    {
                        temp.emplace(fabs(velRat - sbin), sbin);
                    }
                    return temp.begin()->second;
                }
                else
                {
                    return *subcyclingBins.rbegin();
                }

            }
            else
            {
                return *subcyclingBins.begin();
            }
        }

        // return DislocationSegmentIO<dim>::velocityGroup(this->source->get_V(),this->sink->get_V(),maxVelocity);
    }

    template <int dim>
    const typename DislocationSegment<dim>::VectorDim& DislocationSegment<dim>::glidePlaneNormal() const
    {
        return this->glidePlanes().size()==1? (*this->glidePlanes().begin())->unitNormal : zeroVector;
    }
    
    template <int dim>
    typename DislocationSegment<dim>::GlidePlaneContainerType DislocationSegment<dim>::glidePlanes() const
    {
        GlidePlaneContainerType temp;
        for (const auto &ln : this->loopLinks())
        {
            if (ln->periodicPlanePatch())
            {
                temp.emplace(ln->periodicPlanePatch()->glidePlane.get());
            }
        }
        return temp;
    }

    template <int dim>
    typename DislocationSegment<dim>::PlanarMeshFaceContainerType DislocationSegment<dim>::meshFaces() const
    {
        PlanarMeshFaceContainerType temp;
        const PlanarMeshFaceContainerType sourceMeshFaces(this->source->meshFaces());
        const PlanarMeshFaceContainerType sinkMeshFaces(this->sink->meshFaces());

        for (const auto &sourceMeshFace : sourceMeshFaces)
        {
            if (sinkMeshFaces.find(sourceMeshFace) != sinkMeshFaces.end())
            {
                temp.emplace(sourceMeshFace);
            }
        }
        return temp;
    }
    template <int dim>
    bool DislocationSegment<dim>::isOnExternalBoundary() const
    { /*!\returns _isOnExternalBoundarySegment.
       */
        bool _isOnExternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first == face->regionIDs.second)
            {
                _isOnExternalBoundary = true;
            }
        }

        return _isOnExternalBoundary;
    }

    template <int dim>
    bool DislocationSegment<dim>::isOnInternalBoundary() const
    {
        bool _isOnInternalBoundary(false);
        for (const auto &face : meshFaces())
        {
            if (face->regionIDs.first != face->regionIDs.second)
            {
                _isOnInternalBoundary = true;
            }
        }
        return _isOnInternalBoundary;
    }

    template <int dim>
    bool DislocationSegment<dim>::isOnBoundary() const
    {
        return isOnExternalBoundary() || isOnInternalBoundary();
    }
    template <int dim>
    typename DislocationSegment<dim>::GrainContainerType DislocationSegment<dim>::grains() const
    {
        GrainContainerType temp;
        for (const auto &glidePlane : glidePlanes())
        {
            temp.insert(&glidePlane->grain);
        }
        return temp;
    }

    template <int dim>
    typename DislocationSegment<dim>::VectorDim DislocationSegment<dim>::bndNormal() const
    {
        VectorDim _outNormal(VectorDim::Zero());
        for (const auto &face : meshFaces())
        {
            _outNormal += face->outNormal();
        }
        const double _outNormalNorm(_outNormal.norm());
        if (_outNormalNorm > FLT_EPSILON)
        {
            _outNormal /= _outNormalNorm;
        }
        else
        {
            _outNormal.setZero();
        }
        return _outNormal;
    }

    template <int dim>
    std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> DislocationSegment<dim>::parallelAndCoincidentGlidePlanes(const GlidePlaneContainerType &other) const
    {
        std::vector<std::pair<const GlidePlane<dim> *const, const GlidePlane<dim> *const>> pp;

        for (const auto &plane : glidePlanes())
        {
            for (const auto &otherPlane : other)
            {
                //                    if(plane!=otherPlane && plane->n.cross(otherPlane->n).squaredNorm()==0)
                if (plane->n.cross(otherPlane->n).base().squaredNorm() == 0)
                { // parallel planes
                    pp.emplace_back(plane, otherPlane);
                }
            }
        }
        return pp;
    }

    template <int dim>
    typename DislocationSegment<dim>::VectorDim DislocationSegment<dim>::snapToGlidePlanes(const VectorDim &x) const
    {
        GlidePlaneContainerType gps(glidePlanes());
        if(gps.size())
        {
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> N(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,gps.size()));
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,gps.size()));

            int k=0;
            for(const auto& plane : gps)
            {
                N.col(k)=plane->unitNormal;
                P.col(k)=plane->P;
                ++k;
            }

            PlanesIntersection<dim> pInt(N,P,FLT_EPSILON);
            const std::pair<bool,VectorDim> snapped(pInt.snap(x));
            if(snapped.first)
            {
                return snapped.second;
            }
            else
            {
                std::cout<<"glidePlanes.size()="<<gps.size()<<std::endl;
                std::cout<<"N="<<N<<std::endl;
                std::cout<<"P="<<P<<std::endl;
                throw std::runtime_error("DislocationSegment: cannot snap, glidePlanes dont intersect.");
                return snapped.second;
            }
        }
        else
        {
            return x;
        }
    }

    template <int dim>
    typename DislocationSegment<dim>::ConcentrationMatrixType DislocationSegment<dim>::concentrationMatrices(const VectorDim& x,const ClusterDynamicsParameters<dim>& icp) const
    {
        ConcentrationMatrixType M(ConcentrationMatrixType::Zero());
        if(this->grains().size() == 1 && this->chordLength()>FLT_EPSILON)
        {
            const auto& Dinv(icp.invD.at((*this->grains().begin())->grainID));
            const auto& Ddet(icp.detD.at((*this->grains().begin())->grainID));
            const VectorDim bCt(burgers().cross(this->unitDirection()));
            Eigen::Array<double,mSize,1> a(Eigen::Array<double,mSize,1>::Zero());
            Eigen::Array<double,mSize,1> b(Eigen::Array<double,mSize,1>::Zero());
            Eigen::Array<double,mSize,1> c(Eigen::Array<double,mSize,1>::Zero());
            const double bxtSource(bCt.dot(this->source->climbDirection()));
            const double bxtSink  (bCt.dot(this->sink->climbDirection()));
            
            Eigen::Array<double,mSize,1> bias(Eigen::Array<double,mSize,1>::Ones());
            if(this->loopLinks().size()==1)
            {
                const auto loop((*this->loopLinks().begin())->loop);
                if(loop->burgers().dot(loop->rightHandedUnitNormal())>FLT_EPSILON)
                {// vacancy loop
                    bias=icp.discreteDislocationBias.row(0);
                }
                else if(loop->burgers().dot(loop->rightHandedUnitNormal())<-FLT_EPSILON)
                {
                    bias=icp.discreteDislocationBias.row(1);
                }
            }

            for(const auto& shift : this->network().ddBase.periodicShifts)
            {
                for(size_t k=0; k<mSize; k++)
                {
                    a(k) = this->chord().dot( Dinv[k]*this->chord() );
                    b(k) =-2.0*(x+shift-this->source->get_P()).dot( Dinv[k]*this->chord()  );
                    c(k) = (x+shift-this->source->get_P()).dot( Dinv[k]*(x+shift-this->source->get_P()) )+ DislocationFieldBase<dim>::a2/pow(Ddet(k),1.0/3.0);
                }
                const Eigen::Array<double,mSize,1> ba(b/a);
                const Eigen::Array<double,mSize,1> ca(c/a);
                const Eigen::Array<double,mSize,1> sqbca(sqrt(1.0+ba+ca));
                const Eigen::Array<double,mSize,1> sqca(sqrt(ca));
                const Eigen::Array<double,mSize,1> logTerm(log((2.0*sqbca+2.0+ba)/(2.0*sqca+ba)));
                const Eigen::Array<double,mSize,1> I0((1.0+0.5*ba)*logTerm-sqbca+sqca);
                const Eigen::Array<double,mSize,1> I1(     -0.5*ba*logTerm+sqbca-sqca);
                M.col(0)+=bxtSource*this->chordLength()/(4.0*M_PI)*(I0/bias/sqrt(a*Ddet)).matrix();
                M.col(1)+=bxtSink*  this->chordLength()/(4.0*M_PI)*(I1/bias/sqrt(a*Ddet)).matrix();
            }
        }
        return M;
    }

    template <int dim>
    typename DislocationSegment<dim>::ConcentrationVectorType DislocationSegment<dim>::clusterConcentration(const VectorDim& x,const ClusterDynamicsParameters<dim>& icp) const
    {
        ConcentrationVectorType temp(ConcentrationVectorType::Zero());
        const ConcentrationMatrixType cM(concentrationMatrices(x,icp));
        for(int k=0; k<mSize; k++)
        {
            temp(k)=cM.row(k)*(Eigen::Matrix<double,2,1>()<<this->source->climbVelocityScalar(k),this->sink->climbVelocityScalar(k)).finished();
        }
        return temp;
    }

//    template <int dim>
//    const typename DislocationSegment<dim>::MatrixDim DislocationSegment<dim>::I=DislocationSegment<dim>::MatrixDim::Identity();
    
    template <int dim>
    const typename DislocationSegment<dim>::VectorDim DislocationSegment<dim>::zeroVector=DislocationSegment<dim>::VectorDim::Zero();
    
    template <int dim>
    double DislocationSegment<dim>::quadPerLength=0.2;
        
    template <int dim>
    int DislocationSegment<dim>::verboseDislocationSegment=0;
        
    template class DislocationSegment<3>;
}
#endif
