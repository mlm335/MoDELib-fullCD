/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_HEXlattice_cpp_
#define model_HEXlattice_cpp_

#include <HEXlattice.h>
#include <DislocationMobilitySelector.h>
#include <GlidePlaneNoise.h>


namespace model
{

    Eigen::Matrix<double,3,3> HEXlattice<3>::getLatticeBasis()
    {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
      * in units of the crystallographic Burgers vector.
      */
        Eigen::Matrix<double,dim,dim> temp;
        temp << 1.0, 0.5,           0.0,
        /*   */ 0.0, 0.5*sqrt(3.0), 0.0,
        /*   */ 0.0, 0.0,           sqrt(8.0/3.0);
        return temp;
    }

    typename HEXlattice<3>::PlaneNormalContainerType HEXlattice<3>::planeNormals(const PolycrystallineMaterialBase& material,
                                                                                 const Lattice<dim>& lat)
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip plane normals of the HEX lattice
      */
        
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        typedef LatticeVector<dim> LatticeVectorType;
        LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
        LatticeVectorType a3(a2-a1);
        LatticeVectorType  c((VectorDimI()<<0,0,1).finished(),lat);
        
        PlaneNormalContainerType temp;
        if(material.isEnabledPlane("basal"))
        {
            
            const double ISF(TextFileParser(material.materialFile).readScalar<double>("basalISF_SI",true)/(material.mu_SI*material.b_SI));
            const double USF(TextFileParser(material.materialFile).readScalar<double>("basalUSF_SI",true)/(material.mu_SI*material.b_SI));
            const double MSF(TextFileParser(material.materialFile).readScalar<double>("basalMSF_SI",true)/(material.mu_SI*material.b_SI));
            
            const Eigen::Matrix<double,3,2> waveVectors((Eigen::Matrix<double,3,2>()<<0.0, 0.0,
                                                         /*                        */ 0.0, 1.0,
                                                         /*                        */ 1.0,-1.0
                                                         ).finished());
            
            const Eigen::Matrix<double,4,3> f((Eigen::Matrix<double,4,3>()<<0.00,0.0, 0.0,
                                               /*                        */ 0.50,sqrt(3.0)/6.0, ISF,
                                               /*                        */ 0.25,sqrt(3.0)/12.0,USF,
                                               /*                        */ 1.00,sqrt(3.0)/3.0, MSF).finished());
            
            const int rotSymm(3);
            const std::vector<Eigen::Matrix<double,2,1>> mirSymm;
            const Eigen::Matrix<double,2,2> A(GlidePlaneBase(a1,a2,nullptr).localBasis());
            
            std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(A,waveVectors,f,rotSymm,mirSymm));
            temp.emplace(temp.size(),new GlidePlaneBase(a1,a2,gammaSurface));           // basal plane
            
        }
        
        if(material.isEnabledPlane("prismatic"))
        {
            Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors(TextFileParser(material.materialFile).readMatrixCols<double>("prismaticWaveVectors",2,true));
            Eigen::Matrix<double,Eigen::Dynamic,3> fr(TextFileParser(material.materialFile).readMatrixCols<double>("prismaticGammaSurfacePoints",3,true));
            const Eigen::Matrix<double,2,2> A(GlidePlaneBase(a1,c,nullptr).localBasis());
            Eigen::Matrix<double,Eigen::Dynamic,3> f(fr);
            f.block(0,0,fr.rows(),2)=(A*fr.block(0,0,fr.rows(),2).transpose()).transpose();
            f.col(2)=fr.col(2)/(material.mu_SI*material.b_SI);
            const int rotSymm(1);
            std::vector<Eigen::Matrix<double,2,1>> mirSymm;
            mirSymm.push_back((Eigen::Matrix<double,2,1>()<<1.0,0.0).finished()); // symm with respect to local y-axis
            mirSymm.push_back((Eigen::Matrix<double,2,1>()<<0.0,1.0).finished()); // symm with respect to local x-axis
            
            std::shared_ptr<GammaSurface> gammaSurface(new GammaSurface(A,waveVectors,f,rotSymm,mirSymm));
            temp.emplace(temp.size(),new GlidePlaneBase(a1     ,c,gammaSurface));           // prismatic plane
            temp.emplace(temp.size(),new GlidePlaneBase(a3     ,c,gammaSurface));           // prismatic plane
            temp.emplace(temp.size(),new GlidePlaneBase(a2*(-1),c,gammaSurface));           // prismatic plane
        }
        
        if(material.isEnabledPlane("pyramidal"))
        {
            temp.emplace(temp.size(),new GlidePlaneBase(a1     ,a2+c,nullptr));         // pyramidal plane
            temp.emplace(temp.size(),new GlidePlaneBase(a2     ,a3+c,nullptr));         // pyramidal plane
            temp.emplace(temp.size(),new GlidePlaneBase(a3     ,c-a1,nullptr));        // pyramidal plane
            temp.emplace(temp.size(),new GlidePlaneBase(a1*(-1),c-a2,nullptr));       // pyramidal plane
            temp.emplace(temp.size(),new GlidePlaneBase(a2*(-1),c-a3,nullptr));       // pyramidal plane
            temp.emplace(temp.size(),new GlidePlaneBase(a3*(-1),a1+c,nullptr));        // pyramidal plane
        }
        
        return temp;
    }

    typename HEXlattice<3>::SlipSystemContainerType HEXlattice<3>::slipSystems(const PolycrystallineMaterialBase& material,
                                                                               const Lattice<dim>& lat,
                                                                               const PlaneNormalContainerType& plN)
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip systems of the Hexagonal lattice
      */
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
        LatticeVectorType  c((VectorDimI()<<0,0,1).finished(),lat);
        

        DislocationMobilitySelector mobilitySelector;
        std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
        SlipSystemContainerType temp;

        if(material.isEnabledSlipSystem("<a>{basal}"))
        {
            const double dBasal(ReciprocalLatticeDirectionType(a1.cross(a2)).planeSpacing());
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_<a>{basal}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));

            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase.second->planeSpacing()-dBasal)<FLT_EPSILON)
                {
                    const auto& b1(planeBase.second->primitiveVectors.first);
                    const auto& b2(planeBase.second->primitiveVectors.second);
                    const auto b3(b2-b1);
                    
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
        
        if(material.isEnabledSlipSystem("<Shockley>{basal}"))
        {
            const double dBasal(ReciprocalLatticeDirectionType(a1.cross(a2)).planeSpacing());
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_<Shockley>{basal}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));

            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase.second->planeSpacing()-dBasal)<FLT_EPSILON)
                {
                    const auto& b1(planeBase.second->primitiveVectors.first);
                    const auto& b2(planeBase.second->primitiveVectors.second);
                    const auto b3(b2-b1);
                    
                    std::vector<RationalLatticeDirection<3>> slipDirs;
                    slipDirs.emplace_back(Rational(1,3),b1-b3);
                    slipDirs.emplace_back(Rational(1,3),b1-b2);
                    slipDirs.emplace_back(Rational(1,3),b2-b1);
                    slipDirs.emplace_back(Rational(1,3),b2-b3);
                    slipDirs.emplace_back(Rational(1,3),b3-b2);
                    slipDirs.emplace_back(Rational(1,3),b3-b1);

                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                    }
                }
            }
        }
        
        if(material.isEnabledSlipSystem("<a>{prismatic}"))
        {
            const double dPrismatic(ReciprocalLatticeDirectionType(a1.cross(c)).planeSpacing());
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_<a>{prismatic}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));

            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase.second->planeSpacing()-dPrismatic)<FLT_EPSILON)
                {
                    const auto& b(planeBase.second->primitiveVectors.first);

                    
                    std::vector<RationalLatticeDirection<3>> slipDirs;
                    slipDirs.emplace_back(Rational( 1,1),b);
                    slipDirs.emplace_back(Rational(-1,1),b);
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                    }
                }
            }
        }
        
        if(material.isEnabledSlipSystem("<a>{pyramidal}"))
        {
            const double dPyramidal(ReciprocalLatticeDirectionType(a1.cross(a2+c)).planeSpacing());
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_<a>{pyramidal}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));

            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase.second->planeSpacing()-dPyramidal)<FLT_EPSILON)
                {
                    const auto& b(planeBase.second->primitiveVectors.first);

                    
                    std::vector<RationalLatticeDirection<3>> slipDirs;
                    slipDirs.emplace_back(Rational( 1,1),b);
                    slipDirs.emplace_back(Rational(-1,1),b);
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                    }
                }
            }
        }
        
//        const std::string dislocationMobilityTypeBasal(TextFileParser(material.materialFile).readString("dislocationMobilityTypeBasal",true));
//        DislocationMobilitySelector mobilitySelector("HEXbasal");
//        const std::shared_ptr<DislocationMobility> hexMobilityBasal(mobilitySelector.getMobility(dislocationMobilityTypeBasal,material));
//        
//        const std::string dislocationMobilityTypePrismatic(TextFileParser(material.materialFile).readString("dislocationMobilityTypePrismatic",true));
//        DislocationMobilitySelector mobilitySelectorPrismatic("HEXprismatic");
//        const std::shared_ptr<DislocationMobility> hexMobilityPrismatic(mobilitySelectorPrismatic.getMobility(dislocationMobilityTypePrismatic,material));
//        
//        const std::string dislocationMobilityTypePyramidal(TextFileParser(material.materialFile).readString("dislocationMobilityTypePyramidal",true));
//        DislocationMobilitySelector mobilitySelectorPyramidal("HEXpyramidal");
//        const std::shared_ptr<DislocationMobility> hexMobilityPyramidal(mobilitySelectorPyramidal.getMobility(dislocationMobilityTypePyramidal,material));
//        
//        std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
//        
//        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
//        typedef LatticeVector<dim> LatticeVectorType;
//        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
//        LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
//        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
//        LatticeVectorType  c((VectorDimI()<<0,0,1).finished(),lat);
//        
//        const double     dBasal(ReciprocalLatticeDirectionType(a1.cross(a2)).planeSpacing());
//        const double dPrismatic(ReciprocalLatticeDirectionType(a1.cross(c)).planeSpacing());
//        const double dPyramidal(ReciprocalLatticeDirectionType(a1.cross(a2+c)).planeSpacing());
//        
//        SlipSystemContainerType temp;
//        for(const auto& planeBase : plN)
//        {
//            if(std::fabs(planeBase.second->planeSpacing()-dBasal)<FLT_EPSILON)
//            {
//                const auto& a1(planeBase.second->primitiveVectors.first);
//                const auto& a2(planeBase.second->primitiveVectors.second);
//                
//                const auto b1(a1);
//                const auto b2(a2);
//                const auto b3(b2-b1);
//                
//                std::vector<RationalLatticeDirection<3>> slipDirs;
//                
//                if(material.enabledSlipSystems.find("fullBasal")!=material.enabledSlipSystems.end())
//                {
//                    // Full slip systems
//                    slipDirs.emplace_back(Rational( 1,1),b1);
//                    slipDirs.emplace_back(Rational(-1,1),b1);
//                    slipDirs.emplace_back(Rational( 1,1),b2);
//                    slipDirs.emplace_back(Rational(-1,1),b2);
//                    slipDirs.emplace_back(Rational( 1,1),b3);
//                    slipDirs.emplace_back(Rational(-1,1),b3);
//                }
//                
//                if(material.enabledSlipSystems.find("ShockleyBasal")!=material.enabledSlipSystems.end())
//                {
//                    // Shockley partials
//                    slipDirs.emplace_back(Rational(1,3),b1-b3);
//                    slipDirs.emplace_back(Rational(1,3),b1-b2);
//                    slipDirs.emplace_back(Rational(1,3),b2-b1);
//                    slipDirs.emplace_back(Rational(1,3),b2-b3);
//                    slipDirs.emplace_back(Rational(1,3),b3-b2);
//                    slipDirs.emplace_back(Rational(1,3),b3-b1);
//                }
//                
//                for(const auto& slipDir : slipDirs)
//                {
//                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,hexMobilityBasal,planeNoise));
//                }
//                
//            }
//            
//            if(std::fabs(planeBase.second->planeSpacing()-dPrismatic)<FLT_EPSILON)
//            {
//                const auto& b(planeBase.second->primitiveVectors.first);
//                std::vector<RationalLatticeDirection<3>> slipDirs;
//                
//                if(material.enabledSlipSystems.find("fullPrismatic")!=material.enabledSlipSystems.end())
//                {
//                    // Full slip systems
//                    slipDirs.emplace_back(Rational( 1,1),b);
//                    slipDirs.emplace_back(Rational(-1,1),b);
//                }
//                
//                for(const auto& slipDir : slipDirs)
//                {
//                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,hexMobilityPrismatic,planeNoise));
//                }
//            }
//            
//            if(std::fabs(planeBase.second->planeSpacing()-dPyramidal)<FLT_EPSILON)
//            {
//                const auto& b(planeBase.second->primitiveVectors.first);
//                std::vector<RationalLatticeDirection<3>> slipDirs;
//                
//                if(material.enabledSlipSystems.find("fullPyramidal")!=material.enabledSlipSystems.end())
//                {
//                    // Full slip systems
//                    slipDirs.emplace_back(Rational( 1,1),b);
//                    slipDirs.emplace_back(Rational(-1,1),b);
//                }
//                
//                for(const auto& slipDir : slipDirs)
//                {
//                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,hexMobilityPyramidal,planeNoise));
//                }
//            }
//            
//        }
        
        return temp;
    }



    typename HEXlattice<3>::SecondPhaseContainerType HEXlattice<3>::secondPhases(const PolycrystallineMaterialBase& material,
                                                                                 const Lattice<dim>& lat,
                                                                                 const PlaneNormalContainerType& plN)
    {
        SecondPhaseContainerType temp;
        
        for(const std::string& sp : material.enabledSecondPhases)
        {
            throw std::runtime_error("Unnown SecondPhase "+sp+" in HEX crystals.");
        }
        return temp;
    }

} // namespace model
#endif

