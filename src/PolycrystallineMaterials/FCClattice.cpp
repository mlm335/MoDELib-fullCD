/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_FCClattice_cpp_
#define model_FCClattice_cpp_

#include <FCClattice.h>
#include <DislocationMobilitySelector.h>
#include <GlidePlaneNoise.h>

namespace model
{

    Eigen::Matrix<double,3,3> FCClattice<3>::getLatticeBasis()
    {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
      * in units of the crystallographic Burgers vector.
      */
        Eigen::Matrix<double,dim,dim> temp;
        temp << 0.0, 1.0, 1.0,
        /*   */ 1.0, 0.0, 1.0,
        /*   */ 1.0, 1.0, 0.0;
        return temp/sqrt(2.0);
    }

    typename FCClattice<3>::PlaneNormalContainerType FCClattice<3>::planeNormals(const PolycrystallineMaterialBase& material,
                                                                                 const Lattice<dim>& lat)
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip plane normals of the FCC lattice
      */
        LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
        LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),lat);
        
        PlaneNormalContainerType temp;
        
//        const bool enable111planes(true);
//        const bool enable100planes(false);
        
        if(material.isEnabledPlane("{111}"))
        {
            const double ISF(TextFileParser(material.materialFile).readScalar<double>("ISF_SI",true)/(material.mu_SI*material.b_SI));
            const double USF(TextFileParser(material.materialFile).readScalar<double>("USF_SI",true)/(material.mu_SI*material.b_SI));
            const double MSF(TextFileParser(material.materialFile).readScalar<double>("MSF_SI",true)/(material.mu_SI*material.b_SI));
            
            const Eigen::Matrix<double,3,2> waveVectors((Eigen::Matrix<double,3,2>()<<0.0, 0.0,
                                                         /*                        */ 0.0, 1.0,
                                                         /*                        */ 1.0,-1.0).finished());
            
            const Eigen::Matrix<double,4,3> f((Eigen::Matrix<double,4,3>()<<0.00,0.0, 0.0,
                                               /*                        */ 0.50,sqrt(3.0)/6.0, ISF,
                                               /*                        */ 0.25,sqrt(3.0)/12.0,USF,
                                               /*                        */ 1.00,sqrt(3.0)/3.0, MSF).finished());
            
            //        const int rotSymm(3);
            //        const std::vector<Eigen::Matrix<double,2,1>> mirSymm;
            const Eigen::Matrix<double,2,2> A111(GlidePlaneBase(a1,a3,nullptr).localBasis());
            std::shared_ptr<GammaSurface> gammaSurface111(new GammaSurface(A111,waveVectors,f,3,std::vector<Eigen::Matrix<double,2,1>>()));
            
            temp.emplace(temp.size(),new GlidePlaneBase(a1,a3,gammaSurface111));           // is (-1, 1,-1) in cartesian
            temp.emplace(temp.size(),new GlidePlaneBase(a3,a2,gammaSurface111));           // is ( 1,-1,-1) in cartesian
            temp.emplace(temp.size(),new GlidePlaneBase(a2,a1,gammaSurface111));           // is (-1,-1, 1) in cartesian
            temp.emplace(temp.size(),new GlidePlaneBase(a1-a3,a2-a3,gammaSurface111));     // is ( 1, 1, 1) in cartesian
        }
        
        if(   material.isEnabledPlane("{100}")
           || material.isEnabledPlane("{010}")
           || material.isEnabledPlane("{001}"))
        {
            const auto c1(a2+a3-a1);
            const auto c2(a3+a1-a2);
            const auto c3(a1+a2-a3);
            
            // TEMPORARY CONSTANT GSFE FOR 100 planes
            const Eigen::Matrix<double,1,2> waveVectors((Eigen::Matrix<double,1,2>()<<0.0, 0.0).finished());
            const double GSFEval_SI(100.0);
            const Eigen::Matrix<double,1,3> f((Eigen::Matrix<double,1,3>()<<0.5,0.5, GSFEval_SI/(material.mu_SI*material.b_SI)).finished());
            const Eigen::Matrix<double,2,2> A100(GlidePlaneBase(c1,c2,nullptr).localBasis());
            std::shared_ptr<GammaSurface> gammaSurface100(new GammaSurface(A100,waveVectors,f,1,std::vector<Eigen::Matrix<double,2,1>>()));
            
            temp.emplace(temp.size(),new GlidePlaneBase(c1,c2,gammaSurface100));           // is (-1, 1,-1) in cartesian
            temp.emplace(temp.size(),new GlidePlaneBase(c2,c3,gammaSurface100));           // is ( 1,-1,-1) in cartesian
            temp.emplace(temp.size(),new GlidePlaneBase(c3,c1,gammaSurface100));           // is (-1,-1, 1) in cartesian
        }
        return temp;
    }

    typename FCClattice<3>::SlipSystemContainerType FCClattice<3>::slipSystems(const PolycrystallineMaterialBase& material,
                                                                               const Lattice<dim>& lat,
                                                                               const PlaneNormalContainerType& plN)
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip systems of the Hexagonal lattice
      */
        
        DislocationMobilitySelector mobilitySelector;
        SlipSystemContainerType temp;

        if(material.isEnabledSlipSystem("a/2<110>{111}"))
        {// Full dislocations
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/2<110>{111}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
            std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
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
        
        if(material.isEnabledSlipSystem("a/6<112>{111}"))
        {// Shockley dislocations
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/6<112>{111}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
            std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
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
        
        if(material.isEnabledSlipSystem("a/3<112>{111}"))
        {// Kear dislocations
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/6<112>{111}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
            std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
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
                    slipDirs.emplace_back(Rational(2,3),b1-b3);
                    slipDirs.emplace_back(Rational(2,3),b1-b2);
                    slipDirs.emplace_back(Rational(2,3),b2-b1);
                    slipDirs.emplace_back(Rational(2,3),b2-b3);
                    slipDirs.emplace_back(Rational(2,3),b3-b2);
                    slipDirs.emplace_back(Rational(2,3),b3-b1);

                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                    }
                }
            }
        }
        
//        const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("dislocationMobilityType",true));
//        const std::shared_ptr<DislocationMobility> fccMobility(mobilitySelector.getMobility(material,dislocationMobilityType));
//        
//        std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
//        
//        const double d111(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,1.0).finished()).planeSpacing());
//        
//        for(const auto& planeBase : plN)
//        {
//            if(std::fabs(planeBase.second->planeSpacing()-d111)<FLT_EPSILON)
//            {// a {111} plane
//                const auto& a1(planeBase.second->primitiveVectors.first);
//                const auto& a3(planeBase.second->primitiveVectors.second);
//                
//                const auto b1(a1);
//                const auto b2(a3-a1);
//                const auto b3(a3*(-1));
//                
//                std::vector<RationalLatticeDirection<3>> slipDirs;
//                
//                if(   material.enabledSlipSystems.find("Full")!=material.enabledSlipSystems.end()
//                   || material.enabledSlipSystems.find("full")!=material.enabledSlipSystems.end())
//                {
//                    // Full slip systems
//                    slipDirs.emplace_back(Rational( 1,1),b1);
//                    slipDirs.emplace_back(Rational(-1,1),b1);
//                    slipDirs.emplace_back(Rational( 1,1),b2);
//                    slipDirs.emplace_back(Rational(-1,1),b2);
//                    slipDirs.emplace_back(Rational( 1,1),b3);
//                    slipDirs.emplace_back(Rational(-1,1),b3);
//                }
//                if(   material.enabledSlipSystems.find("Shockley")!=material.enabledSlipSystems.end()
//                   || material.enabledSlipSystems.find("shockley")!=material.enabledSlipSystems.end())
//                {
//                    // Shockley partials
//                    slipDirs.emplace_back(Rational(1,3),b1-b3);
//                    slipDirs.emplace_back(Rational(1,3),b1-b2);
//                    slipDirs.emplace_back(Rational(1,3),b2-b1);
//                    slipDirs.emplace_back(Rational(1,3),b2-b3);
//                    slipDirs.emplace_back(Rational(1,3),b3-b2);
//                    slipDirs.emplace_back(Rational(1,3),b3-b1);
//                }
//                if(   material.enabledSlipSystems.find("Kear")!=material.enabledSlipSystems.end()
//                   || material.enabledSlipSystems.find("kear")!=material.enabledSlipSystems.end())
//                {
//                    // Kear partials
//                    slipDirs.emplace_back(Rational(2,3),b1-b3);
//                    slipDirs.emplace_back(Rational(2,3),b1-b2);
//                    slipDirs.emplace_back(Rational(2,3),b2-b1);
//                    slipDirs.emplace_back(Rational(2,3),b2-b3);
//                    slipDirs.emplace_back(Rational(2,3),b3-b2);
//                    slipDirs.emplace_back(Rational(2,3),b3-b1);
//                }
//                
//                for(const auto& slipDir : slipDirs)
//                {
//                    temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,fccMobility,planeNoise));
//                }
//            }
//        }
        return temp;
    }



    typename FCClattice<3>::SecondPhaseContainerType FCClattice<3>::secondPhases(const PolycrystallineMaterialBase& material,
                                                                                 const Lattice<dim>& lat,
                                                                                 const PlaneNormalContainerType& plN)
    {
        
        LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
        LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),lat);
        
        const auto c1(a2+a3-a1);
        const auto c2(a3+a1-a2);
        const auto c3(a1+a2-a3);

        
        SecondPhaseContainerType temp;
        
        for(const std::string& sp : material.enabledSecondPhases)
        {
            if(sp=="L12")
            {
                const double APB(TextFileParser(material.materialFile).readScalar<double>("APB_SI",true)/(material.mu_SI*material.b_SI));
                const double SISF(TextFileParser(material.materialFile).readScalar<double>("SISF_SI",true)/(material.mu_SI*material.b_SI));
                const double CESF(TextFileParser(material.materialFile).readScalar<double>("CESF_SI",true)/(material.mu_SI*material.b_SI));
                const double CISF(TextFileParser(material.materialFile).readScalar<double>("CISF_SI",true)/(material.mu_SI*material.b_SI));
                const double SESF(TextFileParser(material.materialFile).readScalar<double>("SESF_SI",true)/(material.mu_SI*material.b_SI));
                
                const Eigen::Matrix<double,4,2> waveVectors111((Eigen::Matrix<double,4,2>()<<0.0, 0.0,
                                                                /*                        */ 0.0, 1.0,
                                                                /*                        */ 1.0,-1.0,
                                                                //                                                                    /*                        */ 1.0,1.0,
                                                                /*                        */ 2.0, 0.0).finished()); // the factor 0.5 accounts for the fact that A_L12=2*A_fcc
                const Eigen::Matrix<double,6,3> f111((Eigen::Matrix<double,6,3>()<<0.00,0.0, 0.0,
                                                      /*                        */ 0.50,sqrt(3.0)/6.0, CISF,
                                                      /*                        */ 1.00,sqrt(3.0)/3.0,SESF,
                                                      /*                        */ 1.50,sqrt(3.0)/2.0, APB,
                                                      /*                        */ 2.00,2.0*sqrt(3.0)/3.0, SISF,
                                                      /*                        */ 2.50,5.0*sqrt(3.0)/6.0, CESF).finished());
//                const int rotSymm111(3);
//                const std::vector<Eigen::Matrix<double,2,1>> mirSymm111;
                const Eigen::Matrix<double,2,2> A111(2.0*GlidePlaneBase(a1,a3,nullptr).localBasis()); // double the matrix basis
                std::shared_ptr<GammaSurface> gammaSurface111(new GammaSurface(A111,waveVectors111,f111,3,std::vector<Eigen::Matrix<double,2,1>>()));

                
                // TEMPORARY CONSTANT GSFE FOR 100 planes
                const Eigen::Matrix<double,1,2> waveVectors100((Eigen::Matrix<double,1,2>()<<0.0, 0.0).finished());
                const double GSFEval100_SI(50.0);
                const Eigen::Matrix<double,1,3> f100((Eigen::Matrix<double,1,3>()<<0.5,0.5, GSFEval100_SI/(material.mu_SI*material.b_SI)).finished());
                const Eigen::Matrix<double,2,2> A100(2.0*GlidePlaneBase(c1,c2,nullptr).localBasis()); // double the matrix basis
                std::shared_ptr<GammaSurface> gammaSurface100(new GammaSurface(A100,waveVectors100,f100,1,std::vector<Eigen::Matrix<double,2,1>>()));

                
                const double d111(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,1.0).finished()).planeSpacing());
                const double d100(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,0.0,0.0).finished()).planeSpacing());

                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& pn : plN)
                {
                    if(std::abs(pn.second->planeSpacing()-d111)<FLT_EPSILON)
                    {// a 111 plane
                        gsMap.emplace(pn.second.get(),gammaSurface111);
                    }
                    else if(std::abs(pn.second->planeSpacing()-d100)<FLT_EPSILON)
                    {
                        gsMap.emplace(pn.second.get(),gammaSurface100);
                    }
                    else
                    {
                        gsMap.emplace(pn.second.get(),nullptr);
                    }
                }
                
                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("L12",gsMap,plN));
                temp.emplace(sp->sID,sp);
            }
            else
            {
                throw std::runtime_error("Unnown SecondPhase "+sp+" in FCC crystals.");
            }
        }
        return temp;
    }

} // namespace model
#endif

