/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityKinkPair_cpp_
#define _model_DislocationMobilityKinkPair_cpp_

#include <string>
#include <DislocationMobilityKinkPair.h>

namespace model
{

    DislocationMobilityKinkPair::DislocationMobilityKinkPair(const PolycrystallineMaterialBase& material,
                                                             const double& h_in,
                                                             const double& w_in,
                                                             const double& B0_SI,
                                                             const double& B1_SI,
                                                             const double& Bk_SI,
                                                             const double& dH0_eV,
                                                             const double& p_in,
                                                             const double& q_in,
                                                             const double& T0_hom,
                                                             const double& tauC_SI,
                                                             const double& a0_in,
                                                             const double& a1_in,
                                                             const double& a2_in,
                                                             const double& a3_in,
                                                             const double& conjugateAngle_deg):
    /* init */ DislocationMobilityBase("DislocationMobilityKinkPair for "+material.materialName)
    /* init */,kB(this->kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,h(h_in)
    /* init */,w(w_in)
    /* init */,B0(B0_SI*material.cs_SI/(material.mu_SI*material.b_SI))
    /* init */,B1(B1_SI*material.cs_SI/(material.mu_SI*material.b_SI))
    /* init */,Bk(Bk_SI*material.cs_SI/(material.mu_SI*material.b_SI))
    /* init */,dH0(dH0_eV)
    /* init */,p(p_in)
    /* init */,q(q_in)
    /* init */,T0(T0_hom*material.Tm)
    /* init */,tauC(tauC_SI/material.mu_SI)
    /* init */,a0(a0_in)
    /* init */,a1(a1_in)
    /* init */,a2(a2_in)
    /* init */,a3(a3_in)
    /* init */,conjugateAngle(conjugateAngle_deg*std::numbers::pi/180)
    {
    }

    double DislocationMobilityKinkPair::velocity(const MatrixDim& S,
                                                 const VectorDim& b,
                                                 const VectorDim& n,
                                                 const double& T,
                                                 const double& dL,
                                                 const double& dt,
                                                 const std::shared_ptr<StochasticForceGenerator>& sfg) const
    {
        
        const double bNorm=b.norm();
        const VectorDim s = b/bNorm;
        const VectorDim n1 = Eigen::AngleAxisd(conjugateAngle,s)*n;
        
        // Compute components of non-Schmid model
        const double tau=s.transpose()*S*n; // magnitude of resolved shear stress
        const double tauOrt=n.cross(s).transpose()*S*n;
        const double tau1=s.transpose()*S*n1; // resolved schear stress on
        const double tauOrt1=n1.cross(s).transpose()*S*n1;
        const double num(std::fabs(tau+a1*tau1));
        const double den(a0*tauC*sigmoid((a2*tauOrt+a3*tauOrt1)/a0/tauC));
        const double stressRatio(num/den);
        
        // Thermally-activated velocity
        const double dg((stressRatio<1.0)? (std::pow(1.0-std::pow(stressRatio,p),q)-T/T0) : 0.0);
        const double dg1((dg>0.0)? dg : 0.0);
        const double expCoeff(exp(-dH0*dg1/(2.0*kB_eV*T)));
        const double sgm(0.5*sigmoid(-0.5*(0.05-dg1)/0.05));
        const double B(Bk*w/(2.0*h)*(1.0-sgm)+(B0+B1*T)*sgm); //kink-dominated to drag-dominated interpolation
        double vel(std::fabs(tau*bNorm/B*expCoeff));
        
        // Apply stochastic velocity
        if(sfg)
        {
            vel+=sfg->stochasticVelocity(kB,T,B,dL,dt);
        }
        return vel;
    }

    double DislocationMobilityKinkPair::sigmoid(const double & x)
    {
        return 2.0/(1.0+exp(2.0*x));
    }
}
#endif
