/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_GammaSurface_cpp_
#define model_GammaSurface_cpp_

#include <memory>
#include <assert.h>
#include <TerminalColors.h>
#include <LatticeModule.h>
#include <GammaSurface.h>

namespace model
{

GammaSurface::GammaSurface(const MatrixLowerDim& A,
                           const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors,
                           const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
                           const int& rotSymm,
                           const std::vector<Eigen::Matrix<double,lowerDim,1>>& mirSymm) :
/* init */ PeriodicLatticeInterpolant<2>(A,waveVectors,f,rotSymm,mirSymm)
///* init */,G2L(getG2L(a1.cartesian(),a1.cross(a2).cartesian()))
{
    std::cout<<greenColor<<"Creating GammaSurface with local basis\n "<<defaultColor<<this->A<<std::endl;
}

double GammaSurface::misfitEnergy(const VectorLowerDim& b) const
{
    return this->operator()(b);
}


}
#endif
