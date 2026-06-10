/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinClimbSolver_cpp_
#define model_GalerkinClimbSolver_cpp_

#include <deque>


#include <ClusterDynamicsParameters.h>
#include <GalerkinClimbSolver.h>
#include <TerminalColors.h>
#include <EqualIteratorRange.h>
#include <TextFileParser.h>
#include <SparsifiedMatrix.h>

namespace model
{

    template <typename DislocationNetworkType>
    GalerkinClimbSolver<DislocationNetworkType>::GalerkinClimbSolver(const DislocationNetworkType& DN_in,const ClusterDynamics<dim>* const CD_in) :
    /* init */ DislocationClimbSolverBase<DislocationNetworkType>(DN_in,CD_in)
    /* init */,solverType(TextFileParser(CD_in->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("solverType",true))
    /* init */,sparsity_threshold(TextFileParser(CD_in->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("sparsity_threshold",true))
    /* init */,diagonal_compensation(TextFileParser(CD_in->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<double>("diagonal_compensation",true))
    {
        std::cout<<greenBoldColor<<"Creating GalerkinClimbSolver"<<defaultColor<<std::endl;
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::ForceVectorMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterForceKernel(const int& k,const NetworkLinkType& fieldSegment) const
    {
        const Eigen::Matrix<double,1,mSize> deltaC(fieldSegment.quadraturePoint(k).cDD-fieldSegment.quadraturePoint(k).cCD);
        const double u(fieldSegment.quadraturePoint(k).abscissa);
        return  (ForceVectorMatrixType()<<(1.0-u)*deltaC,u*deltaC).finished();
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::ForceVectorMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterForceVector(const NetworkLinkType& fieldSegment) const
    {
        ForceVectorMatrixType F(ForceVectorMatrixType::Zero());
        const auto bxt(fieldSegment.burgers().cross(fieldSegment.chord()));
        NetworkLinkType::QuadratureDynamicType::integrate(fieldSegment.quadraturePoints().size(),this,F,&GalerkinClimbSolver<DislocationNetworkType>::clusterForceKernel,fieldSegment);
        F.row(0)*=bxt.dot(fieldSegment.source->climbDirection());
        F.row(1)*=bxt.dot(fieldSegment.sink->climbDirection());
        return F;
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::StiffnessMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterStiffnessKernel(const int& k,const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const
    {
        StiffnessMatrixType temp(StiffnessMatrixType::Zero());
        if(sourceSegment.grains().size() == 1)
        {
            const auto bxt(fieldSegment.burgers().cross(fieldSegment.chord()));
            const double u(fieldSegment.quadraturePoint(k).abscissa);
            const Eigen::Matrix<double,2,1> m(( Eigen::Matrix<double,2,1>()<<(1.0-u)*bxt.dot(fieldSegment.source->climbDirection()),
                                               /*                               */ u*bxt.dot(fieldSegment.sink->climbDirection())).finished());
            const auto sourceConcentratoinMatrices(sourceSegment.concentrationMatrices(fieldSegment.quadraturePoint(k).r,this->CD->cdp));
            for(int kc=0; kc < mSize; kc++)
            {
                temp.template block<2,2>(kc*2,0)+=m*sourceConcentratoinMatrices.row(kc);
            }
        }
        return temp;
    }

    template <typename DislocationNetworkType>
    typename GalerkinClimbSolver<DislocationNetworkType>::StiffnessMatrixType GalerkinClimbSolver<DislocationNetworkType>::clusterStiffnessMatrix(const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const
    {
        StiffnessMatrixType K(StiffnessMatrixType::Zero());
        NetworkLinkType::QuadratureDynamicType::integrate(fieldSegment.quadraturePoints().size(),this,K,&GalerkinClimbSolver<DislocationNetworkType>::clusterStiffnessKernel,fieldSegment,sourceSegment);
        return K;
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd GalerkinClimbSolver<DislocationNetworkType>::getNodeVelocitiesPipe() const
    {
        Eigen::VectorXd nodeVelocities(Eigen::VectorXd::Zero(this->DN.networkNodes().size()*dim));
        return nodeVelocities;
    }

    template <typename DislocationNetworkType>
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocities()
    {
        std::cout<<", climbSolver "<<std::flush;
        computeClimbScalarVelocitiesBulk();
    }

    template <typename DislocationNetworkType>
    void GalerkinClimbSolver<DislocationNetworkType>::computeClimbScalarVelocitiesBulk()
    {
        std::vector<std::vector<Eigen::Triplet<double>>> lhsT(mSize);
        std::vector<Eigen::VectorXd> Fc(mSize,Eigen::VectorXd::Zero(this->DN.networkNodes().size()));
        std::vector<Eigen::VectorXd> KKc(mSize,Eigen::VectorXd::Zero(this->DN.networkNodes().size()));
        
#ifdef _OPENMP
        const size_t nThreads(omp_get_max_threads());
        std::vector<std::vector<std::vector<Eigen::Triplet<double>>>> lhsTV(nThreads,lhsT);
        std::vector<std::vector<Eigen::VectorXd>> KKcT(nThreads,KKc);
        std::vector<std::vector<Eigen::VectorXd>> FcT(nThreads,Fc);
        const EqualIteratorRange<typename DislocationNetworkType::NetworkLinkContainerType::const_iterator> eir(this->DN.networkLinks().begin(),this->DN.networkLinks().end(),nThreads);
#pragma omp parallel for
        for(size_t thread=0;thread<eir.size();++thread)
        {
            auto& lhsT_ref(lhsTV[thread]);
            auto& Fc_ref(FcT[thread]);
            auto& KKc_ref(KKcT[thread]);
            for(auto fieldLinkIter=eir[thread].first;fieldLinkIter!=eir[thread].second;++fieldLinkIter)
//            for(const auto& fieldLink : this->DN.networkLinks())
            {// sum line-integral part of displacement field per segment
                const auto& fieldLink(*fieldLinkIter);
#else
                auto& lhsT_ref(lhsT);
                auto& Fc_ref(Fc);
                auto& KKc_ref(KKc);
                for(const auto& fieldLink : this->DN.networkLinks())
                {
#endif
////#pragma omp parallel for
//        for(size_t thread=0;thread< eir.size();++thread)
//        {
//            for(auto fieldLinkIter=eir[thread].first;fieldLinkIter!=eir[thread].second;++fieldLinkIter)
////            for(const auto& fieldLink : this->DN.networkLinks())
//            {// sum line-integral part of displacement field per segment
//                const auto& fieldLink(*fieldLinkIter);
                if(   !fieldLink.second.lock()->hasZeroBurgers()
                   && fieldLink.second.lock()->isSessile()
                   && !fieldLink.second.lock()->isBoundarySegment()
                   && !fieldLink.second.lock()->isGrainBoundarySegment()
                   &&  fieldLink.second.lock()->chordLength()>FLT_EPSILON
                   )
                {
                    const size_t i0(fieldLink.second.lock()->source->gID());
                    const size_t i1(fieldLink.second.lock()->  sink->gID());
                    
                    const ForceVectorMatrixType fc(clusterForceVector(*fieldLink.second.lock()));
                    for(int kc=0; kc<mSize; ++kc)
                    {
                        Fc_ref[kc](i0)+=fc(0,kc);
                        Fc_ref[kc](i1)+=fc(1,kc);
                    }
                    

                    // for(const auto& sourceLink : this->DN.networkLinks())
                    // {// sum line-integral part of displacement field per segment
                    //     // if(   !sourceLink.second.lock()->hasZeroBurgers()
                    //     //    && !sourceLink.second.lock()->isBoundarySegment()
                    //     //    && !sourceLink.second.lock()->isGrainBoundarySegment()
                    //     //    &&  sourceLink.second.lock()->chordLength()>FLT_EPSILON
                    //     //    &&  (sourceLink.second.lock()->get_r(0.5) - fieldLink.second.lock()->get_r(0.5)).norm() < sparsity_threshold*(sourceLink.second.lock()->chordLength() + fieldLink.second.lock()->chordLength()) 
                    //     //    )

                    //     // Add nearness condition to skip stiffness contributions from far away segments
                    //     const auto fieldPtr  = fieldLink.second.lock();
                    //     const auto sourcePtr = sourceLink.second.lock();
                    //     const bool sourceIsValid =(!sourcePtr->hasZeroBurgers()
                    //                                 && !sourcePtr->isBoundarySegment()
                    //                                 && !sourcePtr->isGrainBoundarySegment()
                    //                                 &&  sourcePtr->chordLength() > FLT_EPSILON);
                        
                    //     const auto fieldLoopIDs  = fieldPtr->loopIDs();
                    //     const auto sourceLoopIDs = sourcePtr->loopIDs();
                    //     bool sameLoop = false;
                    //     for(const auto& fID : fieldLoopIDs)
                    //     {
                    //         if(sourceLoopIDs.count(fID))
                    //         {
                    //             sameLoop = true;
                    //             break;
                    //         }
                    //     }
                    //     const bool nearEnough = ((sourcePtr->get_r(0.5) - fieldPtr->get_r(0.5)).norm() < sparsity_threshold); //* (sourcePtr->chordLength() + fieldPtr->chordLength()));
                    //     const size_t j0(sourceLink.second.lock()->source->gID());
                    //     const size_t j1(sourceLink.second.lock()->  sink->gID());
                    const auto fieldPtr = fieldLink.second.lock();

                    for(const auto& loopPair : this->DN.loops())
                    {
                        const auto loopPtr = loopPair.second.lock();
                        if(!loopPtr)
                        {
                            continue;
                        }

                        bool sameLoop = false;
                        for(const auto& fID : fieldPtr->loopIDs())
                        {
                            if(fID == loopPtr->sID)
                            {
                                sameLoop = true;
                                break;
                            }
                        }
                        if(sameLoop || loopPtr->distanceTo(*fieldPtr) < sparsity_threshold)
                        {
                            // std::cout << "FIELD SEGMENT (" << i0 << "," << i1 << ") " << "gets LOOP " << loopPtr->sID << " sameLoop=" << sameLoop << std::endl;
                            
                            for(const auto& loopLink : loopPtr->loopLinks())
                            {
                                if(!loopLink->networkLink())
                                {
                                    continue;
                                }

                                const auto sourcePtr = loopLink->networkLink();
                                if(!sourcePtr->hasZeroBurgers() && !sourcePtr->isBoundarySegment() && !sourcePtr->isGrainBoundarySegment() &&  sourcePtr->chordLength() > FLT_EPSILON)
                                {
                                    const size_t j0(sourcePtr->source->gID());
                                    const size_t j1(sourcePtr->sink->gID());
                                    const StiffnessMatrixType kcc(clusterStiffnessMatrix(*fieldPtr,*sourcePtr));
                                    for(int kc=0; kc<mSize; ++kc)
                                    {
                                        const Eigen::Matrix<double,2,2> kccs(kcc.template block<2,2>(2*kc,0));
                                        
                                        if(solverType==0) 
                                        { // lumped solver
                                            KKc_ref[kc](i0)+=0.5*kccs(0,0)+0.5*kccs(0,1);
                                            KKc_ref[kc](j0)+=0.5*kccs(0,0)+0.5*kccs(1,0);
                                            KKc_ref[kc](i1)+=0.5*kccs(1,0)+0.5*kccs(1,1);
                                            KKc_ref[kc](j1)+=0.5*kccs(0,1)+0.5*kccs(1,1);
                                        }
                                        else if(solverType==1 || solverType==2)
                                        { // full solver
                                            const bool forceSym(true);
                                            if (forceSym)
                                            { // Force Symmetry 
                                                lhsT_ref[kc].emplace_back(i0,j0,0.5*kccs(0,0));
                                                lhsT_ref[kc].emplace_back(i0,j1,0.5*kccs(0,1));
                                                lhsT_ref[kc].emplace_back(i1,j0,0.5*kccs(1,0));
                                                lhsT_ref[kc].emplace_back(i1,j1,0.5*kccs(1,1));
                                                lhsT_ref[kc].emplace_back(j0,i0,0.5*kccs(0,0));
                                                lhsT_ref[kc].emplace_back(j1,i0,0.5*kccs(0,1));
                                                lhsT_ref[kc].emplace_back(j0,i1,0.5*kccs(1,0));
                                                lhsT_ref[kc].emplace_back(j1,i1,0.5*kccs(1,1));
                                            }
                                            else
                                            { // Assemble Directly from kcc
                                                lhsT_ref[kc].emplace_back(i0,j0,kccs(0,0));
                                                lhsT_ref[kc].emplace_back(i0,j1,kccs(0,1));
                                                lhsT_ref[kc].emplace_back(i1,j0,kccs(1,0));
                                                lhsT_ref[kc].emplace_back(i1,j1,kccs(1,1));
                                            }
                                        }
                                        else
                                        {
                                            throw std::runtime_error("GalerkinClimbSolver: unknown solverType.");
                                        }
                                    }
                                }
                            }
                        }
                    }
                } 
            }


                        // if(sourceIsValid && (sameLoop || nearEnough))                    
                        // {
                        //     // const size_t j0(sourceLink.second.lock()->source->gID());
                        //     // const size_t j1(sourceLink.second.lock()->  sink->gID());
                        //     const StiffnessMatrixType kcc(clusterStiffnessMatrix(*fieldLink.second.lock(),*sourceLink.second.lock()));
                            
                        //     for(int kc=0; kc<mSize; ++kc)
                        //     {
                        //         const Eigen::Matrix<double,2,2> kccs(kcc.template block<2,2>(2*kc,0));
                                
                        //         if(solverType==0) 
                        //         { // lumped solver
                        //             KKc_ref[kc](i0)+=0.5*kccs(0,0)+0.5*kccs(0,1);
                        //             KKc_ref[kc](j0)+=0.5*kccs(0,0)+0.5*kccs(1,0);
                        //             KKc_ref[kc](i1)+=0.5*kccs(1,0)+0.5*kccs(1,1);
                        //             KKc_ref[kc](j1)+=0.5*kccs(0,1)+0.5*kccs(1,1);
                        //         }
                        //         else if(solverType==1) 
                        //         { // semi-lumped 

                        //             lhsT_ref[kc].emplace_back(i0,j0,0.5*kccs(0,0));
                        //             lhsT_ref[kc].emplace_back(i0,j1,0.5*kccs(0,1));
                        //             lhsT_ref[kc].emplace_back(i1,j0,0.5*kccs(1,0));
                        //             lhsT_ref[kc].emplace_back(i1,j1,0.5*kccs(1,1));
                        //             lhsT_ref[kc].emplace_back(j0,i0,0.5*kccs(0,0));
                        //             lhsT_ref[kc].emplace_back(j1,i0,0.5*kccs(0,1));
                        //             lhsT_ref[kc].emplace_back(j0,i1,0.5*kccs(1,0));
                        //             lhsT_ref[kc].emplace_back(j1,i1,0.5*kccs(1,1));

                        //             // if (i0 <= j0)
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(i0,j0,0.5*kccs(0,0));
                        //             // }
                        //             // if (i0 <= j1)
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(i0,j1,0.5*kccs(0,1));
                        //             // }
                        //             // if (i1 <= j0)
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(i1,j0,0.5*kccs(1,0));
                        //             // }
                        //             // if (i1 <= j1)                                    
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(i1,j1,0.5*kccs(1,1));
                        //             // }
                        //             // if (j0 <= i0)
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(j0,i0,0.5*kccs(0,0));
                        //             // }
                        //             // if (j1 <= i0)
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(j1,i0,0.5*kccs(0,1));
                        //             // }
                        //             // if (j0 <= i1)
                        //             // {   
                        //             //     lhsT_ref[kc].emplace_back(j0,i1,0.5*kccs(1,0));
                        //             // }
                        //             // if (j1 <= i1)
                        //             // {
                        //             //     lhsT_ref[kc].emplace_back(j1,i1,0.5*kccs(1,1));
                        //             // }
                        //         }
                        //         else if(solverType==2)
                        //         { // full solver
                        //             const bool forceSym(true);
                        //             if (forceSym)
                        //             { // Force Symmetry 
                        //                 lhsT_ref[kc].emplace_back(i0,j0,0.5*kccs(0,0));
                        //                 lhsT_ref[kc].emplace_back(i0,j1,0.5*kccs(0,1));
                        //                 lhsT_ref[kc].emplace_back(i1,j0,0.5*kccs(1,0));
                        //                 lhsT_ref[kc].emplace_back(i1,j1,0.5*kccs(1,1));
                        //                 lhsT_ref[kc].emplace_back(j0,i0,0.5*kccs(0,0));
                        //                 lhsT_ref[kc].emplace_back(j1,i0,0.5*kccs(0,1));
                        //                 lhsT_ref[kc].emplace_back(j0,i1,0.5*kccs(1,0));
                        //                 lhsT_ref[kc].emplace_back(j1,i1,0.5*kccs(1,1));
                        //             }
                        //             else
                        //             { // Assemble Directly from kcc
                        //                 lhsT_ref[kc].emplace_back(i0,j0,kccs(0,0));
                        //                 lhsT_ref[kc].emplace_back(i0,j1,kccs(0,1));
                        //                 lhsT_ref[kc].emplace_back(i1,j0,kccs(1,0));
                        //                 lhsT_ref[kc].emplace_back(i1,j1,kccs(1,1));
                        //             }
                        //         }
                        //         else
                        //         {
                        //             throw std::runtime_error("GalerkinClimbSolver: unknown solverType.");
                        //         }
                        //     }
                        // }
                        // else
                        // {
                        //     std::cout
                        //     << "IGNORING interaction : "
                        //     << " field=(" << i0 << "," << i1 << ")"
                        //     << " source=(" << j0 << "," << j1 << ")"
                        //     << " dist="
                        //     << (sourcePtr->get_r(0.5)-fieldPtr->get_r(0.5)).norm()
                        //     << std::endl;
                        // }
        //         }
        //     }
            // }
#ifdef _OPENMP
        }
            for(size_t thread=0;thread< eir.size();++thread)
            {// recombine contributions of different threads
                for(int kc=0; kc<mSize; ++kc)
                {
                    Fc[kc]+=FcT[thread][kc];
                    KKc[kc]+=KKcT[thread][kc];
                    lhsT[kc].insert(lhsT[kc].end(),lhsTV[thread][kc].begin(),lhsTV[thread][kc].end());
                }
            }
#endif
        
        // Eigen::SparseMatrix<double> Kcc(this->DN.networkNodes().size(),this->DN.networkNodes().size());
        std::vector<Eigen::Array<double,1,mSize>> nodeV(this->DN.networkNodes().size(),Eigen::Array<double,1,mSize>::Zero());
        
        for(int kc=0; kc<mSize; ++kc)
        {
            
            if(solverType == 0) 
            { // lumped solver
                for (size_t n=0; n<this->DN.networkNodes().size(); n++)
                {
                    if(std::fabs(KKc[kc](n))>FLT_EPSILON)
                    {
                        nodeV[n](kc)=Fc[kc](n)/KKc[kc](n);
                    }
                }
            }
            else if(solverType == 1) 
            { // semi lumped solver
                
                // Eigen::SparseMatrix<double> Kcc_full(this->DN.networkNodes().size(),this->DN.networkNodes().size());
                // Kcc_full.setFromTriplets(lhsT[kc].begin(), lhsT[kc].end());

                // const int nNodes = static_cast<int>(this->DN.networkNodes().size());
                // const int nnZperRow = 100; 
                // SparsifiedMatrix Kcc(Kcc_full.selfadjointView<Eigen::Upper>(),diagonal_compensation,nnZperRow);

                // std::cout << lhsT[kc].size() << " triplets, "<<Kcc.nonZeros()<<" nonzeros in Kcc matrix for kc="<<kc<<std::endl;
                // std::cout << Kcc.toDense() << std::endl  << std::endl;
                // std::cout << "/nn/n" << std::endl;

                // Eigen::MINRES<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
                // solver.setTolerance(1e-4);
                // solver.setMaxIterations(1e6);
                // solver.compute(Kcc);

                // const Eigen::VectorXd vc = solver.solve(Fc[kc]);
                // for(size_t n=0; n<this->DN.networkNodes().size(); n++)
                // {
                //     nodeV[n](kc)=vc(n);
                // }   

                // // throw std::runtime_error("GalerkinClimbSolver: semi-lumped solver is not implemented yet.");
                Eigen::SparseMatrix<double> Kcc_full(this->DN.networkNodes().size(),this->DN.networkNodes().size());
                Kcc_full.setFromTriplets(lhsT[kc].begin(),lhsT[kc].end());
                
                if(size_t(Kcc_full.rows())!=this->DN.networkNodes().size() || size_t(Kcc_full.cols())!=this->DN.networkNodes().size())
                {
                    throw std::runtime_error("the Stiffness Matrix size is not equal to the node size.");
                }

                // Redue Kcc with threshold 
                Eigen::VectorXd diag = Kcc_full.diagonal();
                std::vector<Eigen::Triplet<double>> reducedTriplets;
                for(int col=0; col<Kcc_full.outerSize(); ++col)
                {
                    for(Eigen::SparseMatrix<double>::InnerIterator it(Kcc_full,col); it; ++it)
                    {
                        const int i = it.row();
                        const int j = it.col();
                        const double Kij = it.value();

                        // if(i > j)
                        // {
                        //     continue; // only upper triangle
                        // }

                        if(i == j)
                        {
                            reducedTriplets.emplace_back(i,j,Kij);
                        }
                        else
                        {
                            // const double threshold = sparsity_threshold * std::sqrt(std::fabs(diag(i)*diag(j)));
                            // const double threshold = sparsity_threshold * std::min(std::fabs(diag(i)),std::fabs(diag(j)));
                            const double threshold = sparsity_threshold * std::fabs(diag(i));
                            if(std::fabs(Kij) < threshold)
                            {
                                reducedTriplets.emplace_back(i,i,diagonal_compensation*std::fabs(Kij));
                                // reducedTriplets.emplace_back(j,j,2*diagonal_compensation*std::fabs(Kij));
                                // std::cout<<"Lumping entry ("<<i<<","<<j<<") with value "<<Kij<<" below threshold "<<threshold<<std::endl;
                            }
                            else
                            {
                                reducedTriplets.emplace_back(i,j,Kij);
                                // reducedTriplets.emplace_back(j,i,Kij);
                            }
                        }
                    }
                }
                // Kcc.setZero();
                Eigen::SparseMatrix<double> Kcc(this->DN.networkNodes().size(),this->DN.networkNodes().size());
                Kcc.setFromTriplets(reducedTriplets.begin(), reducedTriplets.end());

                std::cout << "Lumping Efficiency:" << double( (Kcc_full.nonZeros() - Kcc.nonZeros()) ) / Kcc_full.nonZeros() <<std::endl;

                // Eigen::MINRES<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
                solver.setTolerance(1e-4);
                solver.setMaxIterations(1e6);
                solver.compute(Kcc);

                const Eigen::VectorXd vc = solver.solve(Fc[kc]);
                for(size_t n=0; n<this->DN.networkNodes().size(); n++)
                {
                    nodeV[n](kc)=vc(n);
                }
            }
            else if(solverType == 2) 
            { // full solver
                Eigen::SparseMatrix<double> Kcc(this->DN.networkNodes().size(),this->DN.networkNodes().size());
                Kcc.setFromTriplets(lhsT[kc].begin(), lhsT[kc].end());

                // std::cout << Kcc.toDense() << std::endl;

                std::cout << "EFFICIENCY: " << double(this->DN.networkNodes().size()*this->DN.networkNodes().size() - Kcc.nonZeros()) / (this->DN.networkNodes().size()*this->DN.networkNodes().size()) << std::endl;
                
                if(size_t(Kcc.rows())!=this->DN.networkNodes().size() || size_t(Kcc.cols())!=this->DN.networkNodes().size())
                {
                    throw std::runtime_error("the Stiffness Matrix size is not equal to the node size.");
                }
 
                Eigen::MINRES<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
                solver.setTolerance(1e-4);
                solver.setMaxIterations(1e6);
                solver.compute(Kcc);

                const Eigen::VectorXd vc = solver.solve(Fc[kc]);
                for(size_t n=0; n<this->DN.networkNodes().size(); n++)
                {
                    nodeV[n](kc) = vc(n);
                }     

            }
            else
            {
                throw std::runtime_error("GalerkinClimbSolver: invalid solver type.");
            }
        }

        this->scalarVelocities().resize(this->DN.networkNodes().size(),Eigen::Array<double,1,mSize>::Zero());
        for(size_t k=0; k<this->DN.networkNodes().size();++k)
        {
            this->scalarVelocities()[k]=nodeV[k];
        }
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd GalerkinClimbSolver<DislocationNetworkType>::getNodeVelocitiesBulk() const
    {
        Eigen::VectorXd nodeVelocities(Eigen::VectorXd::Zero(this->DN.networkNodes().size()*dim));
        size_t k=0;
        for(auto& node: this->DN.networkNodes())
        {
            const double vScalarTotal(-1.0*(this->scalarVelocities()[k]*this->CD->cdp.msVector/this->CD->cdp.msVector.abs()).matrix().sum());
            nodeVelocities.template segment<dim>(k*dim)= vScalarTotal*node.second.lock()->climbDirection();
            k++;
        }
        return nodeVelocities;
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd GalerkinClimbSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        return getNodeVelocitiesBulk()+getNodeVelocitiesPipe();
    }

    template class GalerkinClimbSolver<DislocationNetwork<3>>;

}
#endif
