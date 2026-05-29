/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolycrystallineMaterialBase_cpp_
#define model_PolycrystallineMaterialBase_cpp_

#include <string>
#include <TerminalColors.h>
#include <TextFileParser.h>
#include <PolycrystallineMaterialBase.h>

namespace model
{

    PolycrystallineMaterialBase::PolycrystallineMaterialBase(const std::string& fileName,const double& absoluteTemperature) :
    /* init */ materialFile(getMaterialFile(fileName))
    /* init */,materialName(TextFileParser(materialFile).readString("materialName",true))
    /* init */,crystalStructure(TextFileParser(materialFile).readString("crystalStructure",true))
    /* init */,atomsPerUnitCell(TextFileParser(materialFile).readScalar<int>("atomsPerUnitCell",true))
    /* init */,T(absoluteTemperature)
    /* init */,Tm(TextFileParser(materialFile).readScalar<double>("Tm",true))
    /* init */,mu0_SI(TextFileParser(materialFile).readScalar<double>("mu0_SI",true))
    /* init */,mu1_SI(TextFileParser(materialFile).readScalar<double>("mu1_SI",true))
    /* init */,mu_SI(mu0_SI+mu1_SI*T)
    /* init */,nu(TextFileParser(materialFile).readScalar<double>("nu",true))
    /* init */,rho_SI(TextFileParser(materialFile).readScalar<double>("rho_SI",true))
    /* init */,cs_SI(sqrt(mu_SI/rho_SI))
    /* init */,b_SI(TextFileParser(materialFile).readScalar<double>("b_SI",true))
    /* init */,kB(kB_SI/mu_SI/std::pow(b_SI,3))
    /* init */,mu(1.0)
    /* init */,b(1.0)
    /* init */,cs(1.0)
    /* init */,C1(1.0-nu)
    /* init */,C2(1.0/(4.0*std::numbers::pi*C1))
    /* init */,C3(1.0-2.0*nu)
    /* init */,C4(0.5*C2)
    /* init */,enabledSlipSystems(TextFileParser(materialFile).readStringSet("enabledSlipSystems",true))
    /* init */,enabledSecondPhases(TextFileParser(materialFile).readStringSet("enabledSecondPhases",true))
    {
        std::cout<<magentaColor<<"  temperature: T="<<T<<" [K]"<<std::endl;
        std::cout<<magentaColor<<"  units of stress (shear modulus): mu="<<mu_SI<<" [Pa]"<<std::endl;
        std::cout<<magentaColor<<"  units of length (Burgers vector): b="<<b_SI<<" [m]"<<std::endl;
        std::cout<<magentaColor<<"  units of speed (shear-wave speed): cs="<<cs_SI<<" [m/s]"<<std::endl;
        std::cout<<magentaColor<<"  units of time: b/cs="<<b_SI/cs_SI<<" [sec]"<<defaultColor<<std::endl;
        
        std::cout<<"Enabled SlipSystems"<<std::endl;
        for(const auto& str : enabledSlipSystems)
        {
            std::cout<<str<<std::endl;
        }
        
        std::cout<<"Enabled SecondPhases"<<std::endl;
        for(const auto& str : enabledSecondPhases)
        {
            std::cout<<str<<std::endl;
        }
    }

    const std::string& PolycrystallineMaterialBase::getMaterialFile(const std::string& fileName)
    {
        std::cout<<greenBoldColor<<"Reading material file "<<fileName<<defaultColor<<std::endl;
        return fileName;
    }

    bool PolycrystallineMaterialBase::isEnabledPlane(const std::string& planeStr) const
    {
        for(const auto& str : enabledSlipSystems)
        {
            if(str.find(planeStr)!= std::string::npos)
            {
                return true;
            }
        }
        return false;
    }

    bool PolycrystallineMaterialBase::isEnabledSlipSystem(const std::string& ssStr) const
    {
        return enabledSlipSystems.find(ssStr)!=enabledSlipSystems.end();
    }

}
#endif
