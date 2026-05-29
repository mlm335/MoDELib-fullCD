/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilitySelector_cpp_
#define _model_DislocationMobilitySelector_cpp_

#include <numbers>
#include <string>         // std::string
#include <DislocationMobilitySelector.h>
#include <DislocationMobilityEdgeScrew.h>
#include <DislocationMobilityPy.h>
#include <DislocationMobilityInterpolated.h>
#include <DislocationMobilityViscousDrag.h>
#include <DislocationMobilityKinkPair.h>
#include <StrUtilities.h>
#include <TextFileParser.h>
#include <LinearInterpolant.h>
#include <LagrangeInterpolant.h>

namespace model
{

    std::shared_ptr<DislocationMobility> DislocationMobilitySelector::getMobility(const PolycrystallineMaterialBase& material,
                                                                                  const std::string& tag) const
    {
        std::cout<<"tag="<<tag<<std::endl;
        std::stringstream sstream(tag);
        std::string mobilityType;
        sstream >> mobilityType;
        std::cout<<"mobilityType="<<mobilityType<<std::endl;
        
        if(StrUtilities::lowercase(StrUtilities::removeSpaces((mobilityType)))=="none")
        {
            return nullptr;
        }
        else if(StrUtilities::lowercase(StrUtilities::removeSpaces((mobilityType)))=="edgescrew")
        {
            TextFileParser parser(material.materialFile);
            std::string edgeSuffix;
            sstream >> edgeSuffix;
            std::string screwSuffix;
            sstream >> screwSuffix;
            const auto  edgeMobility(getMobilityBase(material,parser.readString("mobility_"+ edgeSuffix)));
            const auto screwMobility(getMobilityBase(material,parser.readString("mobility_"+screwSuffix)));
            return std::shared_ptr<DislocationMobilityEdgeScrew>(new DislocationMobilityEdgeScrew(material,edgeMobility,screwMobility));
        }
        else if(StrUtilities::lowercase(StrUtilities::removeSpaces((mobilityType)))=="python")
        {
            std::string pyPath;
            sstream >> pyPath;
            return std::shared_ptr<DislocationMobilityPy>(new DislocationMobilityPy(material,pyPath));
        }
        else if(StrUtilities::lowercase(StrUtilities::removeSpaces((mobilityType)))=="interpolated")
        {
            TextFileParser parser(material.materialFile);

            std::string interpolationType;
            sstream >> interpolationType;
             ; // periodic extrapolation of period 2pi
            InterpolantBase::MapType vMap;
            std::map<double,std::shared_ptr<DislocationMobilityBase>> mobilityMap;
            std::string angleDeg;
            while (sstream >> angleDeg)
            {
                size_t foundUnderscore(angleDeg.find_last_of('_'));
                if(foundUnderscore != std::string::npos)
                {
                    const double angleRad(std::stod(angleDeg.substr(foundUnderscore+1))*std::numbers::pi/180.0);
                    vMap[angleRad]=0.0;
                    mobilityMap[angleRad]=getMobilityBase(material,parser.readString("mobility_"+ angleDeg));
                }
                else
                {
                    throw std::runtime_error(angleDeg+" NOT in the format label_angleDeg");
                }
                
            }
            
            
            if(StrUtilities::lowercase(StrUtilities::removeSpaces((interpolationType)))=="linear")
            {
                std::shared_ptr<LinearInterpolant> interpolant(new LinearInterpolant(vMap,ExtrapolationMethod(2,std::numbers::pi)));
                return std::shared_ptr<DislocationMobilityInterpolated>(new DislocationMobilityInterpolated(material,interpolant,mobilityMap,false));
            }
            else if(StrUtilities::lowercase(StrUtilities::removeSpaces((interpolationType)))=="logarithmic")
            {
                std::shared_ptr<LinearInterpolant> interpolant(new LinearInterpolant(vMap,ExtrapolationMethod(2,std::numbers::pi)));
                return std::shared_ptr<DislocationMobilityInterpolated>(new DislocationMobilityInterpolated(material,interpolant,mobilityMap,true));
            }
            else if(StrUtilities::lowercase(StrUtilities::removeSpaces((interpolationType)))=="lagrange")
            {
                std::shared_ptr<LagrangeInterpolant> interpolant(new LagrangeInterpolant(vMap,ExtrapolationMethod(2,std::numbers::pi)));
                return std::shared_ptr<DislocationMobilityInterpolated>(new DislocationMobilityInterpolated(material,interpolant,mobilityMap));
            }
            else
            {
                throw std::runtime_error("unknown interpolationType of type "+interpolationType);
                return nullptr;
            }
        }
        else
        {
            throw std::runtime_error("unknown DislocationMobility of type "+mobilityType);
            return nullptr;
        }
    }

    std::shared_ptr<DislocationMobilityBase> DislocationMobilitySelector::getMobilityBase(const PolycrystallineMaterialBase& material,
                                                                                          const std::string& tag) const
    {
        std::stringstream sstream(tag);
        std::string mobilityType;
        sstream >> mobilityType;
        if(StrUtilities::lowercase(StrUtilities::removeSpaces((mobilityType)))=="viscousdrag")
        {
            std::string B0_SI;
            sstream >> B0_SI;
            std::string B1_SI;
            sstream >> B1_SI;
            return std::shared_ptr<DislocationMobilityViscousDrag>(new DislocationMobilityViscousDrag(material,
                                                                                                      std::stod(B0_SI),
                                                                                                      std::stod(B1_SI)));
        }
        else if(StrUtilities::lowercase(StrUtilities::removeSpaces((mobilityType)))=="kinkpair")
        {
            std::string h;
            sstream >> h;
            std::string w;
            sstream >> w;
            std::string B0_SI;
            sstream >> B0_SI;
            std::string B1_SI;
            sstream >> B1_SI;
            std::string Bk_SI;
            sstream >> Bk_SI;
            std::string dH0_eV;
            sstream >> dH0_eV;
            std::string p_in;
            sstream >> p_in;
            std::string q_in;
            sstream >> q_in;
            std::string T0_hom;
            sstream >> T0_hom;
            std::string tauC_SI;
            sstream >> tauC_SI;
            std::string a0_in;
            sstream >> a0_in;
            std::string a1_in;
            sstream >> a1_in;
            std::string a2_in;
            sstream >> a2_in;
            std::string a3_in;
            sstream >> a3_in;
            std::string conjugateAngle_deg;
            sstream >> conjugateAngle_deg;
            
            return std::shared_ptr<DislocationMobilityKinkPair>(new DislocationMobilityKinkPair(material,
                                                                                                std::stod(h),
                                                                                                std::stod(w),
                                                                                                std::stod(B0_SI),
                                                                                                std::stod(B1_SI),
                                                                                                std::stod(Bk_SI),
                                                                                                std::stod(dH0_eV),
                                                                                                std::stod(p_in),
                                                                                                std::stod(q_in),
                                                                                                std::stod(T0_hom),
                                                                                                std::stod(tauC_SI),
                                                                                                std::stod(a0_in),
                                                                                                std::stod(a1_in),
                                                                                                std::stod(a2_in),
                                                                                                std::stod(a3_in),
                                                                                                std::stod(conjugateAngle_deg)));
        }
        else
        {
            throw std::runtime_error("unknown DislocationMobilityBase of type "+mobilityType);
            return nullptr;
        }
    }

}
#endif
