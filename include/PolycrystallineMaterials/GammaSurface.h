/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_GammaSurface_H_
#define model_GammaSurface_H_

#include <memory>
#include <assert.h>
#include <LatticeModule.h>
#include <PeriodicLatticeInterpolant.h>

namespace model
{
    
    struct GammaSurface : public PeriodicLatticeInterpolant<2>
    {
        
        static constexpr int dim=3;
        static constexpr int lowerDim=dim-1;
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,lowerDim,1> VectorLowerDim;
        typedef Eigen::Matrix<size_t,lowerDim,1> VectorLowerDimI;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,lowerDim,lowerDim> MatrixLowerDim;
                
        GammaSurface(const MatrixLowerDim& A,
                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors,
                     const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
                     const int& rotSymm,
                     const std::vector<Eigen::Matrix<double,lowerDim,1>>& mirSymm);
        
        double misfitEnergy(const VectorLowerDim& b) const;

    };
    
}
#endif
