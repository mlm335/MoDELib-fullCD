/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalTraits_H_
#define model_SingleCrystalTraits_H_

#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <map>

#include <TypeTraits.h>
#include <GlidePlaneBase.h>
#include <SlipSystem.h>
#include <SecondPhase.h>

namespace model
{
    template<int dim>
    struct SingleCrystalBase;

    template<int dim>
    struct TypeTraits<SingleCrystalBase<dim>>
    {
        typedef Lattice<dim> LatticeType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::map<size_t,std::shared_ptr<GlidePlaneBase>> PlaneNormalContainerType;
        typedef std::map<size_t,std::shared_ptr<SlipSystem>> SlipSystemContainerType;
        typedef std::map<size_t,std::shared_ptr<SecondPhase<dim>>> SecondPhaseContainerType;
    };
    
}
#endif
