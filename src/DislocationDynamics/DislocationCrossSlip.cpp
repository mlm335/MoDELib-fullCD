/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_cpp_
#define model_DislocationCrossSlip_cpp_

#include <numbers>
#include <memory>

#include <DislocationCrossSlip.h>
#include <CrossSlipModels.h>

namespace model
{
    template <typename DislocationNetworkType>
    std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>> DislocationCrossSlip<DislocationNetworkType>::getModel(const PolycrystallineMaterialBase& material,const DDtraitsIO& traitsIO)
    {
        const int crossSlipModel(TextFileParser(traitsIO.ddFile).readScalar<int>("crossSlipModel",true));
        switch (crossSlipModel)
        {
            case 0:
            {// no cross-slip
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                break;
            }
                
            case 1:
            {// deterministic cross-slip model based on max glide PK force
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(new AthermalCrossSlipModel<DislocationNetworkType>(traitsIO));
                break;
            }
                
            case 2:
            {
                if(material.crystalStructure=="HEX")
                {
                    return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(new EscaigCrossSlipModelHEX<DislocationNetworkType>(material,traitsIO));
                }
                else
                {
                    throw std::runtime_error("Unknown cross-slip model "+std::to_string(crossSlipModel)+" for crystal structure "+material.crystalStructure);
                    return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                }
                break;
            }
                
            default:
            {
                throw std::runtime_error("Unknown cross-slip model "+std::to_string(crossSlipModel));
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                break;
            }
        }
    }

    template <typename DislocationNetworkType>
    DislocationCrossSlip<DislocationNetworkType>::DislocationCrossSlip(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    /* init */,verboseCrossSlip(TextFileParser(DN.ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseCrossSlip",true))
    {
    }

    template <typename DislocationNetworkType>
    void  DislocationCrossSlip<DislocationNetworkType>::findCrossSlipSegments()
    {
        
        if(DN.crossSlipModel)
        {
            std::cout<<"CrossSlip: finding branches "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            csNodes.clear();
            for(const auto& loopIter : DN.loops())
            {
                const auto loop(loopIter.second.lock());
                loop->crossSlipBranches(csNodes);
            }
            std::cout<<csNodes.size()<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
    }

    template <typename DislocationNetworkType>
    void  DislocationCrossSlip<DislocationNetworkType>::execute()
    {
        if(DN.crossSlipModel)
        {
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"CrossSlip: aligning branches"<<std::flush;
//            std::cout<<"csNodes.size()="<<csNodes.size()<<std::endl;
//            std::cout<<"aligning nodes"<<std::endl;
            
            const double dcsMin(20.0);
            
            size_t executed(0);
            CrossSlipBranchContainerType alignedBranches;
            for(const auto& branch : csNodes)
            {// Align all nodes in a brach on the same cross-slip plane
                if(branch.first.size()>=3)
                {
//                    std::cout<<"branch: "<<branch.second<<std::endl;
                    const auto& loop(branch.first.back()->loop());
                    if(loop->glidePlane)
                    {
                        
                        VectorDim midPoint(VectorDim::Zero());
                        for(const auto& node : branch.first)
                        {// Compute the average position of the brach nodes in the superspace
                            midPoint+=node->get_P();
                        }
                        midPoint/=branch.first.size(); // midpoint position in superspace
                        
                        const auto& grain(loop->grain);
                        const auto& crosSlipSystem(grain.slipSystems().at(branch.second));

                        // Compute the heigth of the crossSlip plane closest to midPoint
                        const int dhMin(dcsMin/crosSlipSystem->n.planeSpacing());
                        const double height0(crosSlipSystem->n.closestPlaneIndexOfPoint(midPoint));
                        const long int height(std::round(height0/dhMin)*dhMin);
//                        std::cout<<"dhMin="<<dhMin<<", height0="<<height0<<" -> "<<height<<std::endl;
                        
                        const VectorDim planePoint(height*crosSlipSystem->n.planeSpacing()*crosSlipSystem->unitNormal);
                        PlanePlaneIntersection<dim> ppi(loop->glidePlane->P,
                                                        loop->glidePlane->unitNormal,
                                                        planePoint,
                                                        crosSlipSystem->unitNormal);

                        bool allAligned(true);
                        for(auto& node : branch.first)
                        {// align all nodes in the branch to perfect screw direction
                            const VectorDim oldP(node->networkNode->get_P());
                            const VectorDim deltaP(ppi.P+(node->get_P()-ppi.P).dot(ppi.d)*ppi.d-node->get_P());
                            node->networkNode->trySet_P(node->networkNode->get_P()+deltaP);
                            if((node->networkNode->get_P()-oldP-deltaP).norm()>FLT_EPSILON)
                            {
//                                std::cout<<"node "<<node->tag()<<", cannot move by"<<deltaP.transpose()<<std::endl;
                                allAligned=false;
                            }
                        }
                        if(allAligned)
                        {
                            alignedBranches.push_back(branch);
                        }
                    }
                }
            }
            DN.updateBoundaryNodes(); // After aligning, bnd nodes must be updated
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

            const auto t1= std::chrono::system_clock::now();
            std::cout<<"CrossSlip: creating loops "<<std::flush;
            std::vector<std::shared_ptr<LoopType>> newLoops;
            for(const auto& branch : alignedBranches)
            {
//                if(branch.first.size()>=3)
//                {
//                    std::cout<<"branch: "<<branch.second<<std::endl;
                    const auto& loop(branch.first.back()->loop());
                    if(loop->glidePlane)
                    {
                        std::cout<<"branch:";
                        for(const auto& node : branch.first)
                        {
                            std::cout<<node->tag()<<",";
                        }
                        std::cout<<std::endl;
                        
                        
                        const auto& grain(loop->grain);
                        const auto& crosSlipSystem(grain.slipSystems().at(branch.second)); //THERE COULD BE A PROBLEM HERE WITH SIGN OF B.
//                        std::cout<<"Point="<<branch.first.back()->get_P().transpose()<<std::endl;
//                        std::cout<<"normal="<<crosSlipSystem->unitNormal.transpose()<<std::endl;

                        // Define the cross-slip plane through last NetworkNode in branch
                        LoopNodeType* currentLoopNode(branch.first.back().get());
                        GlidePlaneKey<dim> crossSlipPlaneKey(currentLoopNode->networkNode->get_P(), crosSlipSystem->n);
                        //GlidePlaneKey<dim> crossSlipPlaneKey(branch.first.back()->get_P(), crosSlipSystem->n); TODO: THIS SHOULD BE USING LOOP NODE!
                        const auto crossSlipPlane(DN.ddBase.glidePlaneFactory.getFromKey(crossSlipPlaneKey));
                        auto crossSlipLoop(DN.loops().create(crosSlipSystem->s.cartesian(), crossSlipPlane));
                        const auto periodicCrossSlipPlane(DN.ddBase.periodicGlidePlaneFactory.get(crossSlipPlane->key));
//                        std::cout<<"Constructed crossSlipPlane="<<std::endl;
                        
                        std::deque<std::shared_ptr<LoopNodeType>> crossSlipLoopNodes;
//                        std::cout<<"currentLoopNode="<<currentLoopNode->get_P().transpose()<<std::endl;
//                        std::cout<<"Computing shift"<<std::endl;

//                        VectorDim shift(periodicCrossSlipPlane->findPatch(periodicCrossSlipPlane->referencePlane->localPosition(currentLoopNode->get_P()),VectorDim::Zero()));
                        VectorDim shift(VectorDim::Zero()); // shift of first crossSlipLoopNode is zero, since crossSlipPlane is defined throught this point

//                        std::cout<<"Getting patch"<<std::endl;

                        std::shared_ptr<PeriodicPlanePatch<dim>> crossSlipPatch(periodicCrossSlipPlane->getPatch(shift));
                        
                        // add a new cross-slip node at the end of the branch
//                        std::shared_ptr<NetworkNodeType> newNetNode(DN.networkNodes().create(currentLoopNode->networkNode->get_P()));
                        
//                        bool allContained(true);
//                        std::cout<<"crossSlipPlane->contains("<<newNetNode->tag()<<")? "<<crossSlipPlane->contains(newNetNode->get_P())<<std::endl;
//                        for(auto& node : branch.first)
//                        {
//                            allContained=(allContained && crossSlipPlane->contains(node->networkNode->get_P())); 
//                            std::cout<<"crossSlipPlane->contains("<<node->networkNode->tag()<<")? "<<crossSlipPlane->contains(node->networkNode->get_P())<<std::endl;
//                        }
                        
//                        if(allContained)
//                        {
//                            crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, newNetNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(nullptr,nullptr))); // TODO: THIS SHOULD BE USING LOOP NODE!
                            
                            int crossID(0);
                            while(true)
                            {// kepping adding new cross-slip nodes traversing the branch in reverse order
    //                            std::cout<<"While loop"<<std::endl;
                                
//                                std::cout<<"branchNode "<<currentLoopNode->networkNode->sID<<std::endl;

                                std::shared_ptr<PeriodicPlaneEdge<dim>> edge1(nullptr);
                                std::shared_ptr<PeriodicPlaneEdge<dim>> edge2(nullptr);
                                if(currentLoopNode->periodicPlaneEdge.first)
                                {// a boundary node on one edge
                                    edge1=crossSlipPatch->edgeOnFaces(currentLoopNode->periodicPlaneEdge.first->meshIntersection->faces);
                                    if(!edge1)
                                    {
                                        throw std::runtime_error("crossSlipPatch patch does not have edge1.");
                                    }
                                    if(currentLoopNode->periodicPlaneEdge.second)
                                    {// a boundary node on two edges
                                        edge2=crossSlipPatch->edgeOnFaces(currentLoopNode->periodicPlaneEdge.second->meshIntersection->faces);
                                        if(!edge2)
                                        {
                                            throw std::runtime_error("crossSlipPatch patch does not have edge2.");
                                        }
                                    }
                                }
                                else
                                {// not a boundary node
                                    edge1=nullptr;
                                    edge2=nullptr;
                                }
                                
    //                            std::cout<<"          currentLoopNode->get_P()="<<currentLoopNode->get_P().transpose()<<std::endl;
    //                            std::cout<<"currentLoopNode->networkNode->get_P()="<<currentLoopNode->networkNode->get_P().transpose()<<std::endl;
    //                            std::cout<<"shift="<<shift.transpose()<<std::endl;

    //                            crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
//                                crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));

                                //                            std::cout<<"crossSlipLoopNodes.back()->get_P()="<<crossSlipLoopNodes.back()->get_P().transpose()<<std::endl;

                                if(currentLoopNode==branch.first.back().get())
                                {// place a LoopNode on the cross-slip line
                                    crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
                                }
                                else if(currentLoopNode==branch.first.front().get())
                                {// place LoopNode on the cross-slip line
                                    crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
                                    break;
//                                    // add a new cross-slip node at the beginning of the branch and break
//                                    newNetNode= (currentLoopNode->networkNode==branch.first.back()->networkNode? crossSlipLoopNodes.front()->networkNode : DN.networkNodes().create(currentLoopNode->networkNode->get_P()));
//                                    crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, newNetNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(nullptr,nullptr)));
//                                    break;
                                }
                                else
                                {// place a LoopNode on the cross-slip line, an a new free node
                                    crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, currentLoopNode->networkNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
                                    crossSlipLoopNodes.emplace_front(DN.loopNodes().create(crossSlipLoop, DN.networkNodes().clone(currentLoopNode->networkNode->sID), currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(edge1,edge2)));
                                }
                                
                                if(edge1)
                                {
                                    if(crossID%2 == 0)
                                    {
                                        shift+=edge1->deltaShift;
                                        if(edge2)
                                        {
                                            shift+=edge2->deltaShift;
                                        }
                                        crossSlipPatch=periodicCrossSlipPlane->getPatch(shift);
                                    }
                                    crossID++;
                                }
                                currentLoopNode=currentLoopNode->prev.first;

//                                if(currentLoopNode==branch.first.front().get())
//                                {
//                                    // add a new cross-slip node at the beginning of the branch and break
//                                    newNetNode= (currentLoopNode->networkNode==branch.first.back()->networkNode? crossSlipLoopNodes.front()->networkNode : DN.networkNodes().create(currentLoopNode->networkNode->get_P()));
//                                    crossSlipLoopNodes.emplace_back(DN.loopNodes().create(crossSlipLoop, newNetNode, currentLoopNode->networkNode->get_P()-shift, crossSlipPatch, std::make_pair(nullptr,nullptr)));
//                                    break;
//                                }
//                                else
//                                {
//                                    currentLoopNode=currentLoopNode->prev.first;
//                                }
                            }
                            DN.insertLoop(crossSlipLoop, crossSlipLoopNodes);
                            newLoops.push_back(crossSlipLoop);
                            executed++;
//                        }
//                        else
//                        {
//                            std::cout<<" not all nodes contained in crossSlipPlane "<<std::endl;
//                        }
                    }
//                }
            }
            
//            for(const auto& newLoop : newLoops)
//            {
//                std::cout<<"new loop"<<std::endl;
//                for(const auto& newNodePair : newLoop->nodeSequence())
//                {
//                    std::cout<<newNodePair.first->get_P().transpose()<<","<<newNodePair.first->periodicPlanePatch()->shift.transpose()<<","<<(newNodePair.first->periodicPlaneEdge.first!=nullptr)<<","<<(newNodePair.first->periodicPlaneEdge.second!=nullptr)<<","<<newNodePair.first->networkNode->sID<<std::endl;
//                }
//            }
            
            DN.updateBoundaryNodes(); // After creating loops, bnd nodes must be updated
            
            for(auto& loop :  newLoops)
            {
                for(const auto& link : loop->loopLinks())
                {
                    if(link->networkLink())
                    {
                        link->networkLink()->createQuadraturePoints(false);
                    }
                }
//                loop->computeStackingFaultForces();
                for(const auto& link : loop->loopLinks())
                {
                    if(link->networkLink())
                    {
                        link->networkLink()->updateQuadraturePoints(false);
                    }
                }

            }
            
            std::cout<<executed<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;


//            std::cout<<executed<< "executed "<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
        
        
        
        
    }

    template class DislocationCrossSlip<DislocationNetwork<3>>;
}
#endif
