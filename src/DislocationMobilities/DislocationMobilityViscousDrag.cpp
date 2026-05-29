/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityViscousDrag_cpp_
#define _model_DislocationMobilityViscousDrag_cpp_

#include <cmath>
#include <DislocationMobilityViscousDrag.h>
#include <TextFileParser.h>

namespace model
{
    DislocationMobilityViscousDrag::DislocationMobilityViscousDrag(const PolycrystallineMaterialBase& material,
                                                                   const double& B0_SI,
                                                                   const double& B1_SI) :
    /* init */ DislocationMobilityBase("DislocationMobilityViscousDrag for "+material.materialName)
    /* init */,kB(this->kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,B0(B0_SI*material.cs_SI/(material.mu_SI*material.b_SI))
    /* init */,B1(B1_SI*material.cs_SI/(material.mu_SI*material.b_SI))
    {/*!
      */
    }

//    DislocationMobilityViscousDrag::DislocationMobilityViscousDrag(const std::string& tag,
//                                             const PolycrystallineMaterialBase& material):
//    /* init */ B0(TextFileParser(material.materialFile).readScalar<double>(tag+"B0_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI))
//    /* init */,B1(TextFileParser(material.materialFile).readScalar<double>(tag+"B1_SI",true)*material.cs_SI/(material.mu_SI*material.b_SI))
//    {
//    }


//    std::pair<double,double> DislocationMobilityViscousDrag::velocity(const double& tau,
//                                                           const double& b,
//                                                           const double& T) const
//    {/*!
//      */
//        return velocity(tau*b,T);
//    }
//
//    std::pair<double,double> DislocationMobilityViscousDrag::velocity(const double& tauXb,
//                                                           const double& T) const
//    {/*!
//      */const double B(B0+B1*T);
//        return std::make_pair(std::fabs(tauXb/B),B);
//    }

    double DislocationMobilityViscousDrag::velocity(const MatrixDim& S,
                                         const VectorDim& b,
                                         const VectorDim& n,
                                         const double& T,
                                         const double& dL,
                                         const double& dt,
                                         const std::shared_ptr<StochasticForceGenerator>& sfg) const
    {
        const double B(B0+B1*T);
        const double tauXb(b.transpose()*S*n);
        const double v(std::fabs(tauXb/B));
        if(sfg)
        {
            return v+sfg->stochasticVelocity(kB,T,B,dL,dt);
        }
        else
        {
            return v;
        }
    }

}
#endif


