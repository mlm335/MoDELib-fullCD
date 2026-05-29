
#ifndef model_FirstOrderReaction_H_
#define model_FirstOrderReaction_H_

#include <Eigen/Dense>
#include <EvalFunction.h>
#include <ClusterDynamicsParameters.h>
#include <PolycrystallineMaterialBase.h>

namespace model
{

template <typename MobileTrialFunctionType, typename ImmobileTrialFunctionType>
struct FirstOrderReaction : public EvalFunction<FirstOrderReaction<MobileTrialFunctionType,ImmobileTrialFunctionType>>
{

    constexpr static int rows = MobileTrialFunctionType::rows;
    constexpr static int cols = rows;
    constexpr static int iSize = ImmobileTrialFunctionType::rows;
    constexpr static int mSize = MobileTrialFunctionType::rows;
    constexpr static int dim = MobileTrialFunctionType::dim;

    const ClusterDynamicsParameters<dim>& cdp;
    const PolycrystallineMaterialBase& mat;

    typedef EvalExpression<ImmobileTrialFunctionType> ImmobileEvalFunctionType;
    const ImmobileTrialFunctionType& sink;
    const ImmobileEvalFunctionType sinke;

    const Eigen::Matrix<double,rows,1> inFactors;
    const Eigen::Matrix<double,rows,1> outFactors;

    FirstOrderReaction(const ImmobileTrialFunctionType& sink_in, const ClusterDynamicsParameters<3>& cdp_in, const PolycrystallineMaterialBase& mat_in) :    
    /* init */ cdp(cdp_in),
    /* init */ mat(mat_in),
    /* init */ sink(sink_in),
    /* init */ sinke(sink_in),
    // /* init */ stress(stress_in),
    // /* init */ stresse(eval(stress_in)),
    /* init */ inFactors((1.0/cdp.msVector.abs()).matrix()),
    /* init */ outFactors(cdp.msVector.abs().matrix()) 
    {

    }

    // // /**********************************************************************/
    // Eigen::Array<double,iSize/2,mSize> biasSIPA(const Eigen::Matrix<double,3,3>& ExternalS, const size_t& gID) const
    // {
    //     const Eigen::Matrix<double,3,3> lat(poly.grains().at(gID).latticeBasis);
    //     const Eigen::Matrix<double,3,1> a1(lat.col(0).normalized());
    //     const Eigen::Matrix<double,3,1> a2(lat.col(1).normalized());
    //     const Eigen::Matrix<double,3,1> a3((a2-a1).normalized());
    //     const Eigen::Matrix<double,3,1> c(lat.col(2).normalized());
        
    //     const Eigen::Matrix<double,iSize/2,1> temp((Eigen::Matrix<double,iSize/2,1>()<< -(ExternalS*c).dot(c)*cdp.omega/cdp.kB/cdp.T,  // SIPA c
    //                                                             /*                         */ -(ExternalS*a1).dot(a1)*cdp.omega/cdp.kB/cdp.T,  // SIPA a1
    //                                                             /*                         */ -(ExternalS*a2).dot(a2)*cdp.omega/cdp.kB/cdp.T,  // SIPA a2
    //                                                             /*                         */ -(ExternalS*a3).dot(a3)*cdp.omega/cdp.kB/cdp.T).finished());
        
    //     return exp((temp*cdp.speciesVector.matrix()).array());
    // }
    
    /**********************************************************************/    
    // Eigen::Array<double,iSize/2,mSize> biasDRIFT() const
    // {
    //     const Eigen::Matrix<double,iSize/2,1> temp(Eigen::MatrixXd::Constant(iSize/2,1,1.0));
    //     return (temp*cdp.ZVec.matrix()).array();
    // }
   
    /**********************************************************************/
    Eigen::Array<double,iSize/2,mSize> biasDAD(Eigen::Array<double,1,mSize>& p) const
    {
        Eigen::Array<double,iSize/2,mSize> temp(Eigen::Array<double,iSize/2,mSize>::Zero());
        if(iSize>0)
        {
            for(int i=0; i<iSize/2; i++)
            {
                for(int j=0;j<mSize;j++)
                {
                temp(i,j) = i%(iSize/2)==0 ? p[j] : (p[j]+1.0/p[j]/p[j])/2.0;
                }
            }
        }

        return temp;
    }
    
    /**********************************************************************/
    // Eigen::Matrix<double,iSize,mSize> sinkStrengths(const Eigen::Matrix<double,iSize,1>& sinkValue, const Eigen::Matrix<double,3,3>& ExternalS, const size_t& gID) const
    Eigen::Matrix<double,iSize,mSize> sinkStrengths(const Eigen::Matrix<double,iSize,1>& sinkValue, const size_t& gID) const
    {// Get the stress and the grainID in and compute SIPA in sink strengthes
     // a1, a2, a3 and c can be get in poly.grains.at(grainID).latticeBasis
        
        Eigen::Matrix<double,iSize,1> correctedSinkValue(sinkValue);
        for(int dof=0;dof<iSize/2;dof++)
        { // Corrected sink strength is to adjust for negative values at the barycenter of the quadratic elements where we have a sharp gradient of sink values between two end nodes
            const double cmin = cdp.n_min(dof)*cdp.omega*correctedSinkValue(dof);
            if(correctedSinkValue(iSize/2+dof)<cmin)
            { // defect size cannot be smaller than a minimal value
                correctedSinkValue(iSize/2+dof)=cmin;
            }
        }

        const std::vector<Eigen::Matrix<double,3,3>> Dlocal(cdp.getDlocal());
        Eigen::Array<double,1,mSize> aveD(Eigen::Array<double,1,mSize>::Zero()); // Average diffusion coefficient
        Eigen::Array<double,1,mSize> p(Eigen::Array<double,1,mSize>::Zero());    // Diffusion difference
        for(int k=0;k<mSize;k++)
        {
            aveD(k)=pow(Dlocal[k].determinant(),1.0/3.0);
            p(k)   =pow(Dlocal[k](2,2)/Dlocal[k](0,0),1.0/6.0);
        }
        
        // const Eigen::Array<double,iSize/2,mSize> ZSIPA(biasSIPA(ExternalS,gID));
        // const Eigen::Array<double,iSize/2,mSize> ZDRIFT(biasDRIFT());
        const Eigen::Array<double,iSize/2,mSize> ZDAD(biasDAD(p));
        
        //std::cout<<sinkValue.minCoeff()<<std::endl;
        assert(correctedSinkValue.minCoeff()>=0.0 && "Sink Value is smaller than or equals 0.");
        // std::cout<<"sinkValue: \n"<<sinkValue.transpose()<<std::endl;
        
        // <c> cluster is in pyramidal shape at first
        //Eigen::Array<double,iSize/2,1> loopRadius( (correctedSinkValue.template block<iSize/2,1>(iSize/2,0).array()/correctedSinkValue.template block<iSize/2,1>(0,0).array()/M_PI/cdp.immobileSpeciesBurgersMagnitude.transpose()).sqrt() );
        //Eigen::Matrix<double,iSize/2,1> loopsDensity(2.0*M_PI*loopRadius*correctedSinkValue.template block<iSize/2,1>(0,0).array());
        //loopRadius(0) = cdp.vclusterRadius(correctedSinkValue(iSize/2),correctedSinkValue(0));
        //loopsDensity(0) = 0.2*4.0*M_PI*loopRadius(0)*correctedSinkValue(0);

        const Eigen::Matrix<double,iSize/2,1> disDensity((Eigen::Matrix<double,iSize/2,1>()<<cdp.dislocationSinks(0), //dislocation density c
                                                                                             cdp.dislocationSinks(1), //a1
                                                                                             cdp.dislocationSinks(2), //a2
                                                                                             cdp.dislocationSinks(3)).finished());
        const Eigen::Matrix<double,iSize/2,1> defectDensity(cdp.clusterDensity(correctedSinkValue.template block<iSize/2,1>(iSize/2,0).array(),correctedSinkValue.template block<iSize/2,1>(0,0).array()).transpose().matrix());
        const Eigen::Array<double,iSize/2,mSize> rholD((defectDensity*aveD.matrix()).array());
        const Eigen::Array<double,iSize/2,mSize> rhodD((disDensity*aveD.matrix()).array());
        
        // const Eigen::Array<double,iSize/2,mSize> Zloop = ZSIPA*ZDRIFT*ZDAD;
        const Eigen::Array<double,iSize/2,mSize> Zloop = ZDAD;
        const Eigen::Array<double,iSize/2,mSize> Zpyr(Eigen::Array<double,iSize/2,mSize>::Ones());
        Eigen::Array<double,iSize/2,mSize> Ztotal(Eigen::Array<double,iSize/2,mSize>::Zero());
        
        for(int k=0;k<mSize;k++)
        {
            Ztotal.col(k) = cdp.sigmoidalVectorInterpolation(correctedSinkValue.template block<iSize/2,1>(iSize/2,0).array(),correctedSinkValue.template block<iSize/2,1>(0,0).array(),Zpyr.col(k),Zloop.col(k));
        }
        
        const Eigen::Matrix<double,iSize/2,mSize> loopStrength(Ztotal*rholD); // note that the array multiplication is element-wise
        const Eigen::Matrix<double,iSize/2,mSize> disStrengh(Zloop*rhodD);  // note that the array multiplication is element-wise
        
        Eigen::Matrix<double,iSize,mSize> sStrengh;
        sStrengh << loopStrength,
                    disStrengh;
        
        return sStrengh;
    }
    /**********************************************************************/

    // Eigen::Matrix<double,mSize,mSize> buildR1(const Eigen::Matrix<double,iSize,1>& sinkValue, const Eigen::Matrix<double,3,3>& ExternalS, const size_t& gID) const
    Eigen::Matrix<double,mSize,mSize> buildR1(const Eigen::Matrix<double,iSize,1>& sinkValue, const size_t& gID) const
    {
        Eigen::Matrix<double,iSize,mSize> sStrengh(sinkStrengths(sinkValue,gID));
        
        Eigen::Matrix<double,mSize,mSize>  R1sink(Eigen::Matrix<double,mSize,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            R1sink(k,k)  -= sStrengh.col(k).sum();
        }
        
        return R1sink;
    }
    /**********************************************************************/
    /**********************************************************************/  
    
    
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
    {/*!@param[in] elem the element
      * @param[in] bary the barycentric cooridinate
      *\returns the current stiffness C, which in general is a funciton of C0 and grad(u)
      */
        
        const Eigen::Matrix<double,iSize,1> sinkvalue(sinke(ele,bary));
        const Eigen::Matrix<double,iSize/2,1> loopsDensity((Eigen::Matrix<double,iSize/2,1>()<<2.0*M_PI*sinkvalue(0)*sinkvalue(4),  // Loop density c
                                                                        /*                        */ 2.0*M_PI*sinkvalue(1)*sinkvalue(5),  // Loop density a1
                                                                        /*                        */ 2.0*M_PI*sinkvalue(2)*sinkvalue(6),  // Loop density a2
                                                                        /*                        */ 2.0*M_PI*sinkvalue(3)*sinkvalue(7)).finished());

        // const Eigen::Matrix<double,6,1> tempS(stresse(ele,bary));
        // Eigen::Matrix<double,3,3> ExternalS;
        // ExternalS(0,0)=tempS(0); // s11
        // ExternalS(1,1)=tempS(1); // s22
        // ExternalS(2,2)=tempS(2); // s33
        // ExternalS(1,0)=tempS(3); // s21
        // ExternalS(2,1)=tempS(4); // s32
        // ExternalS(2,0)=tempS(5); // s31
        // ExternalS(0,1)=ExternalS(1,0); //symm
        // ExternalS(1,2)=ExternalS(2,1); //symm
        // ExternalS(0,2)=ExternalS(2,0); //symm


        const size_t gID(*(ele.simplex.regionIDs().begin()));
        // Eigen::Matrix<double,rows,cols> temp(buildR1(sinkvalue,ExternalS,gID));
        Eigen::Matrix<double,rows,cols> temp(buildR1(sinkvalue,gID));
       
        return  outFactors.asDiagonal()*(temp)*inFactors.asDiagonal();
    }

};

} // namespace model

#endif
