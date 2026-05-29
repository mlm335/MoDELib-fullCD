/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobility_h_
#define _model_DislocationMobility_h_

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
 // defines std::cout
#include <TerminalColors.h> // defines mode::cout
#include <PolycrystallineMaterialBase.h>
#include <StaticID.h>
#include <DDtraitsIO.h>
#include <DislocationMobilityBase.h>

namespace model
{
    
    
    
    struct DislocationMobility : public StaticID<DislocationMobility>
    {
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        const std::string name;

        
        DislocationMobility(const std::string& name_in) ;
        
        virtual ~DislocationMobility(){};
        
        virtual double velocity(const MatrixDim& S,
                                const VectorDim& b,
                                const VectorDim& , // xi
                                const VectorDim& n,
                                const double& T,
                                const double& dL,
                                const double& dt,
                                const std::shared_ptr<StochasticForceGenerator>& sfg) const =0 ;
        
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T) const ;
        
    };
    
}
#endif
