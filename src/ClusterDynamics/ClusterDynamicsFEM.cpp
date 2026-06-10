/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_cpp_
#define model_ClusterDynamicsFEM_cpp_

#include <cmath>
#include <ClusterDynamicsFEM.h>
#include <ExternalAndInternalBoundary.h>
#include <Fix.h>

namespace model
{

    template<int dim>
    FluxMatrix<dim>::FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */cdp(cdp_in)
    {
        
    }

    template<int dim>
    const typename FluxMatrix<dim>::MatrixType FluxMatrix<dim>::operator() (const ElementType& elem, const BaryType&) const
    {

        const size_t grainID(elem.simplex.region->regionID);
        MatrixType D(MatrixType::Zero());
        for(size_t k=0;k<mSize;++k)
        {
            D.template block<dim,dim>(k*dim,k*dim)=-cdp.D.at(grainID)[k]/cdp.omega;
        }
        
        // for(size_t k=0;k<mSize;++k)
        // {
        //     std::cout
        //         << " species " << k
        //         << " D trace = "
        //         << cdp.D.at(grainID)[k].trace()
        //         << std::endl;
        // }

        return D;
    }

    template struct FluxMatrix<3>;

template<int dim>
InvDscaling<dim>::InvDscaling(const ClusterDynamicsParameters<dim>& cdp_in) :
/* init */cdp(cdp_in)
{
    
}

template<int dim>
const typename InvDscaling<dim>::MatrixType InvDscaling<dim>::operator() (const ElementType& elem, const BaryType&) const
{
    const size_t grainID(elem.simplex.region->regionID);
    MatrixType InvDscaling(MatrixType::Zero());
    for(size_t k=0;k<mSize;++k)
    {
        InvDscaling(k,k)=3.0/(cdp.D.at(grainID)[k].trace()/cdp.omega);
//        InvDscaling(k,k)=(cdp.D.at(grainID)[k].trace()/cdp.omega)/3.0;
//        InvDscaling(k,k)=1.0;

    }
    return InvDscaling;
}

template struct InvDscaling<3>;


    template<int dim>
    ClusterDynamicsFEM<dim>::ClusterDynamicsFEM(DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */ ddBase(ddBase_in)
    /* init */,cdp(cdp_in)
    /* init */,iDs(cdp)
    /* init */,mobileClusters(ddBase.fe->template trial<'m',mSize>())
    /* init */,mobileGrad(grad(mobileClusters))
    /* init */,mobileFlux(FluxMatrix<dim>(cdp)*mobileGrad)
    /* init */,immobileClusters(ddBase.fe->template trial<'i',iSize>())
    /* init */,nodeListInternalExternal(ddBase.isPeriodicDomain ? -1 : ddBase.fe->template createNodeList<ExternalAndInternalBoundary>())
    /* init */,mobileClustersIncrement(ddBase.fe->template trial<'d',mSize>())
    /* init */,dV(ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
//    /* init */,mBWF((test(this->mobileGrad),-ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,mBWF((test(grad(iDs*mobileClusters)),-ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,dmBWF((test(grad(iDs*mobileClustersIncrement)),-ddBase.poly.Omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)   
    /* init */,mSolver(true,FLT_EPSILON)
    /* init */,solverInitialized(false)
    /* init */,cascadeGlobalProduction(((test(iDs*this->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())
    // /* init */,cascadeGlobalProduction(((test(this->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())

    {
        mobileClustersIncrement.setConstant(Eigen::Matrix<double,mSize,1>::Zero());
        mobileClusters.setConstant(cdp.equilibriumMobileConcentration(0.0).matrix().transpose());
        // immobileClusters.setConstant(Eigen::Matrix<double,iSize,1>::Zero());
        immobileClusters.setConstant(cdp.initLoopSinks);    
        immobileClusterRate = Eigen::VectorXd::Zero(immobileClusters.gSize());    
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveMobileClusters(const bool hasDiscreteLoops)
    {
        std::cout<<", mobile solver, "<<std::flush;
        mobileClusters=mSolver.solve(cascadeGlobalProduction);
        const bool useImmobileClusters(!hasDiscreteLoops && iSize > 0);

        if(this->cdp.computeReactions)
        {
            const double cTol(1e-5);
            double cError(1.0);
            while(cError>cTol)
            {
                const auto R1((this->cdp.R1cd).eval());
                auto bWF_R1((test(iDs*mobileClustersIncrement),R1*(-1.0*mobileClustersIncrement))*dV); // THIS SHOULD BE STORED SINCE IT IS ALWAYS THE SAME
                auto lWF_R1((test(iDs*mobileClustersIncrement),eval(R1*mobileClusters))*dV);

                FirstOrderReaction<MobileTrialType,ImmobileTrialType> R1sink(immobileClusters,this->cdp,ddBase.poly);
                auto bWF_R1sink((test(iDs*mobileClustersIncrement),R1sink*(-1.0*mobileClustersIncrement))*dV); 
                auto lWF_R1sink((test(iDs*mobileClustersIncrement),eval(R1sink*mobileClusters))*dV);
                
                auto lWF_R1sink2((test(mobileClustersIncrement),eval(R1sink*mobileClusters))*dV);

                SecondOrderReaction<MobileTrialType> R2(mobileClusters,this->cdp);
                auto bWF_R2((test(iDs*mobileClustersIncrement),R2*(-1.0*mobileClustersIncrement))*dV);
                auto lWF_R2((test(iDs*mobileClustersIncrement),eval(R2*(0.5*mobileClusters)))*dV);

                Eigen::SparseMatrix<double,Eigen::RowMajor> AcIR;
                AcIR.resize(mobileClustersIncrement.gSize(),mobileClustersIncrement.gSize());

                std::vector<Eigen::Triplet<double>> globalTripletsR(useImmobileClusters? (bWF_R1+bWF_R2+bWF_R1sink).globalTriplets() : (bWF_R1+bWF_R2).globalTriplets());
                AcIR.setFromTriplets(globalTripletsR.begin(),globalTripletsR.end());
                
                MobileReactionSolverType rSolver(false,FLT_EPSILON);

                if(useImmobileClusters)
                {
                    rSolver.compute(dmBWF+bWF_R1+bWF_R2+bWF_R1sink); 
                }
                else
                {
                    rSolver.compute(dmBWF+bWF_R1+bWF_R2);
                }
                mobileClustersIncrement=rSolver.solve(cascadeGlobalProduction-mSolver.getA()*mobileClusters.dofVector() + (useImmobileClusters? (lWF_R1+lWF_R2+lWF_R1sink).globalVector() : (lWF_R1+lWF_R2).globalVector()) );
                
                Eigen::MatrixXd cOld(mobileClusters.dofVector());
                cOld.resize(mSize,mobileClusters.gSize()/mSize);
                mobileClusters += mobileClustersIncrement.dofVector();
                
                Eigen::MatrixXd cNew(mobileClusters.dofVector());
                cNew.resize(mSize,mobileClusters.gSize()/mSize);
                
                const Eigen::VectorXd absErr((cNew-cOld).rowwise().norm());
                const Eigen::VectorXd cNewNorm((cNew.rowwise().norm().array()+1.e-50).matrix());
                const Eigen::VectorXd relErr((absErr.array()/cNewNorm.array()).matrix());
                
                cError=relErr.maxCoeff();//aError/cInorm;
                if(false)
                {
                    std::cout<<"max values="<<cNew.rowwise().maxCoeff().transpose()<<std::endl;
                    std::cout<<"min values="<<cNew.rowwise().minCoeff().transpose()<<std::endl;
                    std::cout<<"absolute errors="<<absErr.transpose()<<std::endl;
                    std::cout<<"solution norms="<<cNewNorm.transpose()<<std::endl;
                    std::cout<<"relative error="<<relErr.transpose()<<std::endl;
                }
                std::cout<<"convergenceError="<<cError<<std::endl;
            }
        }
        // Find immobile rate // did this below
        // Find diffusive-displacement rate // did this below, i think
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveImmobileClusters()
    {
        if (iSize > 0)
        {
            std::cout<<", immobile solver, "<<std::flush;
            ImmobileSinkRate<MobileTrialType,ImmobileTrialType> immobileRate(mobileClusters,immobileClusters,this->cdp,ddBase.poly);
            auto lWFsink((test(immobileClusters),immobileRate)*dV);
            SpatialODESolver iSolver(immobileClusters,dV,false,1e-4);
            immobileClusterRate = iSolver.solve(lWFsink.globalVector());
            std::cout<<"convergenceError="<<iSolver.error()<<std::endl;
        }
    }


    template<int dim>
    void ClusterDynamicsFEM<dim>::updateImmobileClusters(const double dt)
    {
        TrialBase<ImmobileTrialType>::dofVector()+=immobileClusterRate * dt;
        for(size_t nodeID=0; nodeID<immobileClusters.nodeSize(); nodeID++)
        {
            for(int dof=0;dof<iSize/2;dof++)
            { // Corrected sink strength is to adjust for negative values at the barycenter of the quadratic elements where we have a sharp gradient of sink values between two end nodes
                // rmin = r_min * omega * N
                const double cmin = this->cdp.n_min(dof)*this->cdp.omega*TrialBase<ImmobileTrialType>::dofVector()(iSize*nodeID+dof);
                if(TrialBase<ImmobileTrialType>::dofVector()(iSize*nodeID+iSize/2+dof)<cmin) //  Nc, Na1, Na2, Na3, c, ca1, ca2, ca3
                { // Loop size cannot be smaller than a minimal value
                    TrialBase<ImmobileTrialType>::dofVector()(iSize*nodeID+iSize/2+dof)=cmin;                
                }
            }
        }
    }


    template<int dim>
    void ClusterDynamicsFEM<dim>::solve(const bool hasDiscreteLoops)
    {
        solveMobileClusters(hasDiscreteLoops);
        const bool useImmobileClusters(!hasDiscreteLoops && iSize > 0);
        if(useImmobileClusters)
        {
            solveImmobileClusters();
        }
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeSolver()
    {
        std::cout<<" decomposing"<<std::flush;
        std::array<bool,mSize> allComps;
        for(int k=0;k<mSize;++k)
        {
            allComps[k]=true;
        }
        
        mobileClusters.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
        mobileClustersIncrement.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
        mSolver.compute(mBWF); // call this after assigning the BCs
        solverInitialized=true;
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {    
        if(size_t(configIO.cdMatrix().size())==mobileClusters.gSize()+immobileClusters.gSize())
        {
            const size_t nNodes(mobileClusters.fe().nodes().size());
            mobileClusters=configIO.cdMatrix().block(0,0,nNodes,mSize).transpose().reshaped(mobileClusters.gSize(),1);
            immobileClusters=configIO.cdMatrix().block(0,mSize,nNodes,iSize).transpose().reshaped(immobileClusters.gSize(),1);

        }
        else
        {
            if(configIO.cdMatrix().size())
            {
                throw std::runtime_error("ClusterDynamics: TrialFunctions initializatoin size mismatch");
            }
        }        
    }


    template<int dim>
    typename ClusterDynamicsFEM<dim>::VectorDim ClusterDynamicsFEM<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        Eigen::Matrix<double,dim*mSize,1> speciesFlux(Eigen::Matrix<double,dim*mSize,1>::Zero());
        if(node)
        {
            speciesFlux=eval(this->mobileFlux)(*node);
        }
        else
        {
            if(ele)
            {
                speciesFlux=eval(this->mobileFlux)(*ele,ele->simplex.pos2bary(x));
            }
            else
            {
                speciesFlux=eval(this->mobileFlux)(x,guess);
            }
        }
        Eigen::Matrix<double,dim,1> netFlux(Eigen::Matrix<double,dim,1>::Zero());
        for(int i=0; i<mSize; i++)
        {
            const int mSgn(this->cdp.msVector(i)/std::abs(this->cdp.msVector(i)));
            netFlux += speciesFlux.template block<dim,1>(i*dim,0)*mSgn;
        }
        return netFlux;
    }


    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::betaP(const ElementType& ele, const BaryType& bary) const
    {
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector;
        const Eigen::Matrix<double,iSize,1> sinkValue(eval(immobileClusters)(ele,bary));
        Eigen::Matrix<double,iSize,1> correctedSinkValue(sinkValue);
        for(int dof=0;dof<iSize/2;dof++)
        {
            const double cmin = cdp.n_min(dof)*cdp.omega*correctedSinkValue(dof);
            if(correctedSinkValue(iSize/2+dof)<cmin)
            { // Loop size cannot be smaller than a minimal value
                correctedSinkValue(iSize/2+dof)=cmin;
            }
        }
        const size_t gID(ele.simplex.region->regionID);
        const Eigen::Matrix<double,dim,dim> lat(ddBase.poly.grains.at(gID)->latticeBasis);
        const Eigen::Matrix<double,dim,iSize/2> localBurgers(lat*cdp.immobileSpeciesBurgers);
        const Eigen::Array<double,1,iSize/2> cr(cdp.clusterRadius(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0)));
        const Eigen::Matrix<double,1,dim> localx(lat.col(0).normalized());
        const Eigen::Matrix<double,1,dim> localy((2.0*lat.col(1)-lat.col(0)).normalized());
        const Eigen::Matrix<double,1,dim> localz(lat.col(2).normalized());
        for(int k=0; k<iSize/2; k++)
        { // N A b ^b*^b
            Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
            const double Area=M_PI*std::pow(cr(k),2);
            const Eigen::Matrix<double,dim,dim> highValue(cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*Area*localBurgers.col(k)*localBurgers.col(k).normalized().transpose());
            const Eigen::Matrix<double,dim,dim> lowValue(cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*cdp.delVPyramid*(localx.transpose()*localx+localy.transpose()*localy+localz.transpose()*localz)/3.0);
            tempStrain+=cdp.sigmoidalMatrixInterpolation(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0),lowValue,highValue,k).matrix();
            strainVector.push_back(tempStrain);
        }
        return strainVector;
    }


    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::betaV(const ElementType& ele, const BaryType& bary) const
    {
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector;
        const Eigen::Matrix<double,iSize,1> sinkValue(eval(immobileClusters)(ele,bary));
        Eigen::Matrix<double,iSize,1> correctedSinkValue(sinkValue);
        for(int dof=0;dof<iSize/2;dof++)
        {
        const double cmin = cdp.n_min(dof)*cdp.omega*correctedSinkValue(dof);
            if(correctedSinkValue(iSize/2+dof)<cmin)
            { // Loop size cannot be smaller than a minimal valu
                correctedSinkValue(iSize/2+dof)=cmin;
            }
        }
        const size_t gID(ele.simplex.region->regionID);
        const Eigen::Matrix<double,dim,dim> lat(ddBase.poly.grains.at(gID)->latticeBasis);
        const Eigen::Matrix<double,dim,iSize/2> localBurgers(lat*cdp.immobileSpeciesBurgers);
        const Eigen::Matrix<double,1,dim> localx(lat.col(0).normalized());
        const Eigen::Matrix<double,1,dim> localy((2.0*lat.col(1)-lat.col(0)).normalized());
        const Eigen::Matrix<double,1,dim> localz(lat.col(2).normalized());
        const Eigen::Array<double,1,iSize/2> cr(cdp.clusterRadius(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0)));
        for(int k=0; k<iSize/2; k++)
        { // N A b ^b*^b
            Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
            const double Area=M_PI*std::pow(cr(k),2);
            const Eigen::Matrix<double,dim,dim> highValue(cdp.immobileSpeciesRelRelaxVol(k)*cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*Area*localBurgers.col(k)*localBurgers.col(k).normalized().transpose());
            const Eigen::Matrix<double,dim,dim> lowValue(cdp.immobileSpeciesRelRelaxVol(k)*cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*cdp.delVPyramid*(localx.transpose()*localx+localy.transpose()*localy+localz.transpose()*localz)/3.0);
            tempStrain+=cdp.sigmoidalMatrixInterpolation(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0),lowValue,highValue,k).matrix();
            strainVector.push_back(tempStrain);
        }
       // Add up the contribution from the mobile species
        Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
        const Eigen::Matrix<double,mSize,1> inFactors((1.0/cdp.msVector.abs()).matrix());
        const Eigen::Matrix<double,mSize,1> localConcentration(inFactors.asDiagonal()*eval(mobileClusters)(ele,bary)); // Local concentration
        for(int k=0; k<mSize; k++)
        {
            tempStrain += 1.0/3.0*cdp.msVector(k)*localConcentration(k)*cdp.msRelRelaxVol(k)*Eigen::Matrix<double,dim,dim>::Identity();
        }
        strainVector.push_back(tempStrain);
        return strainVector;
    }
    

    /* plastic distortion */
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaPKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const
    {
        const auto bary(BarycentricTraits<dim>::x2l(a1));

        const std::vector<Eigen::Matrix<double,dim,dim>> bps(betaP(ele,bary));
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());  
        for(const Eigen::Matrix<double,dim,dim>& bp : bps)
        {
            temp += bp*ele.absJ(bary);
        }
        return temp; 
        // return betaP(ele,bary)[0]*ele.absJ(bary);
    }
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaP() const
    {
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());
        // dV.integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaPKernel);    // integrate through the volume
        for(size_t k=0; k<dV.size(); ++k)
        {
            VolumeIntegrationDomainType::QuadratureType::integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaPKernel,dV.element(k));
        }
        return temp/ddBase.poly.mesh.volume();
    }

    /* volumetric distortion */   
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaVKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const
    {
        const auto bary(BarycentricTraits<dim>::x2l(a1));

        const std::vector<Eigen::Matrix<double,dim,dim>> bvs(betaV(ele,bary));
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());  
        for(const Eigen::Matrix<double,dim,dim>& bv : bvs)
        {
            temp += bv*ele.absJ(bary);
        }
        return temp; 
        // return betaP(ele,bary)[0]*ele.absJ(bary);
    }
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaV() const
    {
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());
        // dV.integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaVKernel);    // integrate through the volume
        for(size_t k=0; k<dV.size(); ++k)
        {
            VolumeIntegrationDomainType::QuadratureType::integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaVKernel,dV.element(k));
        }
        return temp/ddBase.poly.mesh.volume();
    }


    template struct ClusterDynamicsFEM<3>;

}
#endif



//template<int dim>
//void ClusterDynamicsFEM<dim>::applyBoundaryConditions()
//{
//
//    const auto& nodesInternalExternal(mobileClusters.fe().nodeList(nodeListInternalExternal));
//
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for(size_t k=0;k<nodesInternalExternal.size();++k)
//    {
//        const auto& node(nodesInternalExternal[k]);
//        const auto outNormal(node->outNormal()); // used to compute traction
//        const MatrixDim sigma(microstructures.stress(node->P0,node,nullptr,nullptr));
//        const double normalTraction(outNormal.dot(sigma*outNormal));
//        const auto bndConcentration(this->cdp.boundaryMobileConcentration(sigma.trace(),normalTraction));
//
//        VectorMSize otherConcentration(VectorMSize::Zero());
//        for(const auto& microstructure : this->microstructures)
//        {
//            if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
//            {// not the ClusterDynamics physics
//                otherConcentration += microstructure->mobileConcentration(node->P0,node,nullptr,nullptr);
//            }
//        }
//        for(int k=0;k<mSize;++k)
//        {
//            mobileClusters.dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
//        }
//    }
//}
