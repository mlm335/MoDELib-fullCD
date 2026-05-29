/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_BCClattice_cpp_
#define model_BCClattice_cpp_

#include <BCClattice.h>
#include <DislocationMobilitySelector.h>
#include <GlidePlaneNoise.h>

namespace model
{

    Eigen::Matrix<double,3,3> BCClattice<3>::getLatticeBasis()
    {/*!\returns The matrix of lattice vectors (cartesian cooridinates in columns),
      * in units of the crystallographic Burgers vector.
      */
        Eigen::Matrix<double,dim,dim> temp;
        temp << -1.0,  1.0,  1.0,
        /*   */  1.0, -1.0,  1.0,
        /*   */  1.0,  1.0, -1.0;
        return temp/sqrt(3.0);
    }

    typename BCClattice<3>::PlaneNormalContainerType BCClattice<3>::planeNormals(const PolycrystallineMaterialBase& material,
                                                                                 const Lattice<dim>& lat)
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip plane normals of the BCC lattice
      */
        
        LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
        LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),lat);
        LatticeVectorType  y((VectorDimI()<<-1,-1,-1).finished(),lat);
        
        PlaneNormalContainerType temp;
        
        if(   material.isEnabledPlane("{110}")
           || material.isEnabledPlane("{101}")
           || material.isEnabledPlane("{011}"))
        {
            temp.emplace(temp.size(),new GlidePlaneBase(a3,a1,nullptr));
            temp.emplace(temp.size(),new GlidePlaneBase( y,a2,nullptr));
            temp.emplace(temp.size(),new GlidePlaneBase(a2,a3,nullptr));
            temp.emplace(temp.size(),new GlidePlaneBase( y,a1,nullptr));
            temp.emplace(temp.size(),new GlidePlaneBase(a1,a2,nullptr));
            temp.emplace(temp.size(),new GlidePlaneBase( y,a3,nullptr));
        }
        return temp;
    }

    typename BCClattice<3>::SlipSystemContainerType BCClattice<3>::slipSystems(const PolycrystallineMaterialBase& material,
                                                                               const Lattice<dim>& lat,
                                                                               const PlaneNormalContainerType& plN)
    {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
      * the slip systems of the Hexagonal lattice
      */
        
        DislocationMobilitySelector mobilitySelector;
        
        SlipSystemContainerType temp;
        
        if(material.isEnabledSlipSystem("a/2<111>{110}"))
        {
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a/2<111>{110}",true));
            const std::shared_ptr<DislocationMobility> mobility(mobilitySelector.getMobility(material,dislocationMobilityType));
            std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
            const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
            
            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase.second->planeSpacing()-d110)<FLT_EPSILON)
                {// a {110} plane
                    const auto& a1(planeBase.second->primitiveVectors.first);
                    const auto& a3(planeBase.second->primitiveVectors.second);
                    
                    std::vector<RationalLatticeDirection<3>> slipDirs;
                    slipDirs.emplace_back(Rational( 1,1),a1);
                    slipDirs.emplace_back(Rational(-1,1),a1);
                    slipDirs.emplace_back(Rational( 1,1),a3);
                    slipDirs.emplace_back(Rational(-1,1),a3);
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility,planeNoise));
                    }
                }
            }
        }
        
        if(material.isEnabledSlipSystem("a<100>{110}"))
        {
            const std::string dislocationMobilityType(TextFileParser(material.materialFile).readString("mobility_a<100>{110}",true));
            const std::shared_ptr<DislocationMobility> mobility110(mobilitySelector.getMobility(material,dislocationMobilityType));
            std::shared_ptr<GlidePlaneNoise> planeNoise(new GlidePlaneNoise(material));
            const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
            
            for(const auto& planeBase : plN)
            {
                if(std::fabs(planeBase.second->planeSpacing()-d110)<FLT_EPSILON)
                {// a {110} plane
                    const auto& a1(planeBase.second->primitiveVectors.first);
                    const auto& a3(planeBase.second->primitiveVectors.second);
                    
                    std::vector<RationalLatticeDirection<3>> slipDirs;
                    slipDirs.emplace_back(Rational( 1,1),a1+a3);
                    slipDirs.emplace_back(Rational(-1,1),a1+a3);
                    
                    for(const auto& slipDir : slipDirs)
                    {
                        temp.emplace(temp.size(),new SlipSystem(*planeBase.second, slipDir,mobility110,planeNoise));
                    }
                }
            }
        }
        
        
        return temp;
    }

    typename BCClattice<3>::SecondPhaseContainerType BCClattice<3>::secondPhases(const PolycrystallineMaterialBase& material,
                                                                                 const Lattice<dim>& lat,
                                                                                 const PlaneNormalContainerType& plN)
    {       LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
        LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),lat);
        
        SecondPhaseContainerType temp;
        
        for(const std::string& sp : material.enabledSecondPhases)
        {
            if(sp=="chi")
            {
                Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors110(TextFileParser(material.materialFile).readMatrixCols<double>("chiWaveVectors110",2,true));
                Eigen::Matrix<double,Eigen::Dynamic,3> fr110(TextFileParser(material.materialFile).readMatrixCols<double>("chiGammaSurfacePoints110",3,true));
                const int rotSymm110(1);
                const Eigen::Matrix<double,2,2> A110(GlidePlaneBase(a3,a1,nullptr).localBasis()*(Eigen::Matrix<double,2,2>()<<2, 1, 0, 1).finished());
                std::vector<Eigen::Matrix<double,2,1>> mirSymm110 = {
                    (A110.col(1)).normalized(),
                    (A110.col(0)-A110.col(1)).normalized()
                };
                
                Eigen::Matrix<double,Eigen::Dynamic,3> f110(fr110);
                f110.block(0,0,fr110.rows(),2) = (A110 * fr110.block(0,0,fr110.rows(),2).transpose()).transpose();
                f110.col(2) = fr110.col(2) / (material.mu_SI * material.b_SI);
                
                std::shared_ptr<GammaSurface> gammaSurface110(new GammaSurface(A110,waveVectors110,f110,rotSymm110,mirSymm110));
                const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& pn : plN)
                {
                    if(std::abs(pn.second->planeSpacing()-d110)<FLT_EPSILON)
                    {
                        gsMap.emplace(pn.second.get(),gammaSurface110);
                    }
                }
                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("chi",gsMap,plN));
                temp.emplace(sp->sID,sp);
            }
            else if(sp=="sigma")
            {
                Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors110(TextFileParser(material.materialFile).readMatrixCols<double>("sigmaWaveVectors110",2,true));
                Eigen::Matrix<double,Eigen::Dynamic,3> fr110(TextFileParser(material.materialFile).readMatrixCols<double>("sigmaGammaSurfacePoints110",3,true));
                const int rotSymm110(1);
                const Eigen::Matrix<double,2,2> A110(GlidePlaneBase(a3,a1,nullptr).localBasis()*(Eigen::Matrix<double,2,2>()<<2, 1, 0, 1).finished());
                std::vector<Eigen::Matrix<double,2,1>> mirSymm110 = {
                    (A110.col(1)).normalized(),
                    (A110.col(0)-A110.col(1)).normalized()
                };
                
                Eigen::Matrix<double,Eigen::Dynamic,3> f110(fr110);
                f110.block(0,0,fr110.rows(),2) = (A110 * fr110.block(0,0,fr110.rows(),2).transpose()).transpose();
                f110.col(2) = fr110.col(2) / (material.mu_SI * material.b_SI);
                
                std::shared_ptr<GammaSurface> gammaSurface110(new GammaSurface(A110,waveVectors110,f110,rotSymm110,mirSymm110));
                const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
                for(const auto& pn : plN)
                {
                    if(std::abs(pn.second->planeSpacing()-d110)<FLT_EPSILON)
                    {
                        gsMap.emplace(pn.second.get(),gammaSurface110);
                    }
                    else
                    {
                        gsMap.emplace(pn.second.get(),nullptr);
                    }
                }
                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("sigma",gsMap,plN));
                temp.emplace(sp->sID,sp);
            }
            
            else
            {
                throw std::runtime_error("Unknown SecondPhase "+sp+" in BCC crystals.");
            }
        }
        return temp;
    }

//    typename BCClattice<3>::SecondPhaseContainerType BCClattice<3>::secondPhases(const PolycrystallineMaterialBase& material,
//                                                                                 const Lattice<dim>& lat,
//                                                                                 const PlaneNormalContainerType& plN)
//    {       LatticeVectorType a1((VectorDimI()<<1,0,0).finished(),lat);
//        LatticeVectorType a2((VectorDimI()<<0,1,0).finished(),lat);
//        LatticeVectorType a3((VectorDimI()<<0,0,1).finished(),lat);
//
//        SecondPhaseContainerType temp;
//
//        for(const std::string& sp : material.enabledSecondPhases)
//        {
//            if(sp=="chi")
//            {
//                Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors110(TextFileParser(material.materialFile).readMatrixCols<double>("chiWaveVectors110",2,true));
//                Eigen::Matrix<double,Eigen::Dynamic,3> fr110(TextFileParser(material.materialFile).readMatrixCols<double>("chiGammaSurfacePoints110",3,true));
//                const int rotSymm110(1);
//                const Eigen::Matrix<double,2,2> A110(2.0*GlidePlaneBase(a3,a1,nullptr).localBasis());
//                std::vector<Eigen::Matrix<double,2,1>> mirSymm110 = {
//                    (A110.col(0)+A110.col(1)).normalized(),
//                    (A110.col(0)-A110.col(1)).normalized()
//                };
//
//                Eigen::Matrix<double,Eigen::Dynamic,3> f110(fr110);
//                f110.block(0,0,fr110.rows(),2) = (A110 * fr110.block(0,0,fr110.rows(),2).transpose()).transpose();
//                f110.col(2) = fr110.col(2) / (material.mu_SI * material.b_SI);
//
//                std::shared_ptr<GammaSurface> gammaSurface110(new GammaSurface(A110,waveVectors110,f110,rotSymm110,mirSymm110));
//                const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
//                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
//                for(const auto& pn : plN)
//                {
//                    if(std::abs(pn.second->planeSpacing()-d110)<FLT_EPSILON)
//                    {
//                        gsMap.emplace(pn.second.get(),gammaSurface110);
//                    }
//                }
//                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("chi",gsMap,plN));
//                temp.emplace(sp->sID,sp);
//            }
//            else if(sp=="sigma")
//            {
//                Eigen::Matrix<double,Eigen::Dynamic,2> waveVectors110(TextFileParser(material.materialFile).readMatrixCols<double>("sigmaWaveVectors110",2,true));
//                Eigen::Matrix<double,Eigen::Dynamic,3> fr110(TextFileParser(material.materialFile).readMatrixCols<double>("sigmaGammaSurfacePoints110",3,true));
//                const int rotSymm110(1);
//                const Eigen::Matrix<double,2,2> A110(2.0*GlidePlaneBase(a3,a1,nullptr).localBasis());
//                std::vector<Eigen::Matrix<double,2,1>> mirSymm110 = {
//                    (A110.col(0)+A110.col(1)).normalized(),
//                    (A110.col(0)-A110.col(1)).normalized()
//                };
//
//                Eigen::Matrix<double,Eigen::Dynamic,3> f110(fr110);
//                f110.block(0,0,fr110.rows(),2) = (A110 * fr110.block(0,0,fr110.rows(),2).transpose()).transpose();
//                f110.col(2) = fr110.col(2) / (material.mu_SI * material.b_SI);
//
//                std::shared_ptr<GammaSurface> gammaSurface110(new GammaSurface(A110,waveVectors110,f110,rotSymm110,mirSymm110));
//                const double d110(lat.reciprocalLatticeDirection(lat.C2G*(VectorDimD()<<1.0,1.0,0.0).finished()).planeSpacing());
//                std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
//                for(const auto& pn : plN)
//                {
//                    if(std::abs(pn.second->planeSpacing()-d110)<FLT_EPSILON)
//                    {
//                        gsMap.emplace(pn.second.get(),gammaSurface110);
//                    }
//                }
//                std::shared_ptr<SecondPhase<dim>> sp(new SecondPhase<dim>("sigma",gsMap,plN));
//                temp.emplace(sp->sID,sp);
//            }
//
//            else
//            {
//                throw std::runtime_error("Unnown SecondPhase "+sp+" in BCC crystals.");
//            }
//        }
//        return temp;
//    }

} // namespace model
#endif
