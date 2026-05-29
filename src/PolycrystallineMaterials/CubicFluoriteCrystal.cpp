/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CubicFluoriteCrystal_cpp_
#define model_CubicFluoriteCrystal_cpp_

#include <CubicFluoriteCrystal.h>
#include <DislocationMobilitySelector.h>
#include <GlidePlaneNoise.h>


namespace model
{

Eigen::Matrix<double,3,3> CubicFluoriteCrystal<3>::getLatticeBasis()
{/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
  * in units of the crystallographic Burgers vector.
  */
    Eigen::Matrix<double,dim,dim> temp;
    temp << 0.0, 1.0, 1.0,
    /*   */ 1.0, 0.0, 1.0,
    /*   */ 1.0, 1.0, 0.0;
    return temp/sqrt(2.0);
}

typename CubicFluoriteCrystal<3>::PlaneNormalContainerType CubicFluoriteCrystal<3>::planeNormals(const PolycrystallineMaterialBase& material,
                                                                                                 const Lattice<dim>& lat)
{/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
  * the slip plane normals of the FCC lattice
  */
    LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
    LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
    LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),lat);
    LatticeVectorType x1((VectorDimI()<<-1,1,1).finished(),lat);
    LatticeVectorType x2((VectorDimI()<<1,-1,1).finished(),lat);
    LatticeVectorType x3((VectorDimI()<<1,1,-1).finished(),lat);
    
    PlaneNormalContainerType temp;
    
    if(   material.isEnabledPlane("{111}"))
    {// 111 planes
        temp.emplace(temp.size(),new GlidePlaneBase(a1,a3,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a3,a2,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a2,a1,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a1-a3,a2-a3,nullptr));
    }
    
    if(   material.isEnabledPlane("{011}")
       || material.isEnabledPlane("{101}")
       || material.isEnabledPlane("{110}"))
    {// 110 planes
        temp.emplace(temp.size(),new GlidePlaneBase(a3,x3,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a2-a1,x3,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a2,x2,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a1-a3,x2,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a1,x1,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a3-a2,x1,nullptr));
    }
    
    if(   material.isEnabledPlane("{100}")
       || material.isEnabledPlane("{010}")
       || material.isEnabledPlane("{001}"))
    {// 100 planes
        temp.emplace(temp.size(),new GlidePlaneBase(a1,a2-a3,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a2,a3-a1,nullptr));
        temp.emplace(temp.size(),new GlidePlaneBase(a3,a1-a2,nullptr));
    }
    
    return temp;
}

typename CubicFluoriteCrystal<3>::SlipSystemContainerType CubicFluoriteCrystal<3>::slipSystems(const PolycrystallineMaterialBase& material,
                                                                                               const Lattice<dim>& lat,
                                                                                               const PlaneNormalContainerType& plN)
{/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
  * the slip systems of the Hexagonal lattice
  */
    
//    const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("dislocationMobilityType",true));
//    DislocationMobilitySelector mobilitySelector("BCC");
//    const std::shared_ptr<DislocationMobilityBase> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
    
    
    DislocationMobilitySelector mobilitySelector;
    std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
    SlipSystemContainerType temp;
        
    if(material.isEnabledSlipSystem("a/2<110>{111}"))
    {
        
        const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/2<110>{111}",true));
        const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
        const double d111(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,1.0).finished()).planeSpacing());
        for(const auto& planeBase : plN)
        {
            if(std::fabs(planeBase.second->planeSpacing()-d111)<FLT_EPSILON)
            {// a {111} plane
                const auto& a1(planeBase.second->primitiveVectors.first);
                const auto& a3(planeBase.second->primitiveVectors.second);
                
                const auto b1(a1);
                const auto b2(a3-a1);
                const auto b3(a3*(-1));
                
                std::vector<RationalLatticeDirection<3>> slipDirs;
                slipDirs.emplace_back(Rational( 1,1),b1);
                slipDirs.emplace_back(Rational(-1,1),b1);
                slipDirs.emplace_back(Rational( 1,1),b2);
                slipDirs.emplace_back(Rational(-1,1),b2);
                slipDirs.emplace_back(Rational( 1,1),b3);
                slipDirs.emplace_back(Rational(-1,1),b3);
                
                for(const auto& slipDir : slipDirs)
                {
                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                }
            }
        }
    }
    
    
    if(material.isEnabledSlipSystem("a/2<110>{110}"))
    {
        const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/2<110>{110}",true));
        const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
        const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
        for(const auto& planeBase : plN)
        {
            if(std::fabs(planeBase.second->planeSpacing()-d110)<FLT_EPSILON)
            {// a {111} plane
                const auto b1(planeBase.second->primitiveVectors.first);
                
                std::vector<RationalLatticeDirection<3>> slipDirs;
                slipDirs.emplace_back(Rational( 1,1),b1);
                slipDirs.emplace_back(Rational(-1,1),b1);
                
                for(const auto& slipDir : slipDirs)
                {
                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                }
            }
        }
    }
    
    if(material.isEnabledSlipSystem("a/2<110>{100}"))
    {
        const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/2<110>{100}",true));
        const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
        const double d100(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,0.0,0.0).finished()).planeSpacing());
        for(const auto& planeBase : plN)
        {
            if(std::fabs(planeBase.second->planeSpacing()-d100)<FLT_EPSILON)
            {// a {111} plane
                const auto& b1(planeBase.second->primitiveVectors.first);
                const auto& b2(planeBase.second->primitiveVectors.second);
                
                std::vector<RationalLatticeDirection<3>> slipDirs;
                slipDirs.emplace_back(Rational( 1,1),b1);
                slipDirs.emplace_back(Rational(-1,1),b1);
                slipDirs.emplace_back(Rational( 1,1),b2);
                slipDirs.emplace_back(Rational(-1,1),b2);
                
                for(const auto& slipDir : slipDirs)
                {
                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                }
            }
        }
    }
    
//    for(const auto& planeBase : plN)
//    {
//        if(std::fabs(planeBase.second->planeSpacing()-d111)<FLT_EPSILON)
//        {// a {111} plane
//            const auto& a1(planeBase.second->primitiveVectors.first);
//            const auto& a3(planeBase.second->primitiveVectors.second);
//            
//            const auto b1(a1);
//            const auto b2(a3-a1);
//            const auto b3(a3*(-1));
//            
//            std::vector<RationalLatticeDirection<3>> slipDirs;
//            
//            if(   material.enabledSlipSystems.find("full<110>{111}")!=material.enabledSlipSystems.end())
//            {
//                // Full slip systems
//                slipDirs.emplace_back(Rational( 1,1),b1);
//                slipDirs.emplace_back(Rational(-1,1),b1);
//                slipDirs.emplace_back(Rational( 1,1),b2);
//                slipDirs.emplace_back(Rational(-1,1),b2);
//                slipDirs.emplace_back(Rational( 1,1),b3);
//                slipDirs.emplace_back(Rational(-1,1),b3);
//            }
//            
//            for(const auto& slipDir : slipDirs)
//            {
//                temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
//            }
//        }
//        
//        if(std::fabs(planeBase.second->planeSpacing()-d110)<FLT_EPSILON)
//        {// a {110} plane
//            const auto& a1(planeBase.second->primitiveVectors.first);
//            const auto b1(a1);
//            
//            std::vector<RationalLatticeDirection<3>> slipDirs;
//            
//            if(   material.enabledSlipSystems.find("full<110>{110}")!=material.enabledSlipSystems.end())
//            {
//                // Full slip systems
//                slipDirs.emplace_back(Rational( 1,1),b1);
//                slipDirs.emplace_back(Rational(-1,1),b1);
//            }
//            
//            for(const auto& slipDir : slipDirs)
//            {
//                temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
//            }
//        }
//        
//        if(std::fabs(planeBase.second->planeSpacing()-d100)<FLT_EPSILON)
//        {// a {100} plane
//            const auto& a1(planeBase.second->primitiveVectors.first);
//            const auto& a2(planeBase.second->primitiveVectors.second);
//            
//            const auto b1(a1);
//            const auto b2(a2);
//            
//            
//            std::vector<RationalLatticeDirection<3>> slipDirs;
//            
//            if(   material.enabledSlipSystems.find("full<110>{100}")!=material.enabledSlipSystems.end())
//            {
//                // Full slip systems
//                slipDirs.emplace_back(Rational( 1,1),b1);
//                slipDirs.emplace_back(Rational(-1,1),b1);
//                slipDirs.emplace_back(Rational( 1,1),b2);
//                slipDirs.emplace_back(Rational(-1,1),b2);
//            }
//            
//            for(const auto& slipDir : slipDirs)
//            {
//                temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
//            }
//        }
//        
//    }
    return temp;
}



typename CubicFluoriteCrystal<3>::SecondPhaseContainerType CubicFluoriteCrystal<3>::secondPhases(const PolycrystallineMaterialBase& material,
                                                                                                 const Lattice<dim>& lat,
                                                                                                 const PlaneNormalContainerType& plN)
{
    
    
    SecondPhaseContainerType temp;
    
    
    return temp;
}

} // namespace model
#endif

