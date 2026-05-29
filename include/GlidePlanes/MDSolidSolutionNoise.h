/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDSolidSolutionNoise_h
#define model_MDSolidSolutionNoise_h

#include <cmath>
#include <random>

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>

#include <PolycrystallineMaterialBase.h>
#include <NoiseTraits.h>
#include <GlidePlaneNoiseBase.h>

namespace model
{

    struct MDSolidSolutionNoise:  public GlidePlaneNoiseBase<2>
    {
        typedef typename NoiseTraits<2>::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraits<2>::COMPLEX COMPLEX;
        typedef typename NoiseTraits<2>::GridSizeType GridSizeType;
        typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
        typedef typename NoiseTraits<2>::NoiseType NoiseType;
        typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;

        COMPLEX *Rk_xz;
        COMPLEX *Rk_yz;
        const REAL_SCALAR a_cai;

        // variable energy correction factor from zero padding
        REAL_SCALAR eNormFactor;
        // fftScaleFactor is the Scaling factor for FFT
        REAL_SCALAR fftScaleFactor;

        MDSolidSolutionNoise(const PolycrystallineMaterialBase& mat,
                                    const std::string& tag,
                                    const std::string& correlationFile_xz, const std::string& correlationFile_yz,
                                    const int& seed, const GridSizeType& gridSize, const GridSpacingType& gridSpacing,
                                    const Eigen::Matrix<double,2,2>& latticeBasis,
                                    const double& a_Cai_in);
        
        
        std::array<COMPLEX,2> kCorrelations(const Eigen::Matrix<double,3,1>& kv,const Eigen::Matrix<int,3,1>& kvID) const override;
        static void SolidSolutionCorrelationReader(const std::string& correlationFile, REAL_SCALAR* Rr, const int& NR);
        static GridSizeType readVTKfileDimension(const char* fname);
    };

}
#endif
