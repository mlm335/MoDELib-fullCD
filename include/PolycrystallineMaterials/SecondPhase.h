/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SecondPhase_H_
#define model_SecondPhase_H_

#include <memory>
#include <map>
#include <assert.h>

#include <LatticeModule.h>
#include <GlidePlaneBase.h>
#include <SlipSystem.h>
#include <StaticID.h>

namespace model
{
    template<int dim>
    struct SecondPhase : public StaticID<SecondPhase<dim>>
    {
        typedef const Eigen::Matrix<double,dim,1> VectorDim;
        typedef std::map<size_t,std::shared_ptr<GlidePlaneBase>> PlaneNormalContainerType;

        const std::string name;
        const std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>> gsMap;
        const PlaneNormalContainerType& planes;
                
        SecondPhase(const std::string&,
                    const std::map<const GlidePlaneBase*,std::shared_ptr<GammaSurface>>&,
                    const PlaneNormalContainerType&);
        double misfitEnergy(const VectorDim& s,const GlidePlaneBase* const gpb) const;
        double misfitEnergy(const VectorDim& s,const size_t&) const;
        std::shared_ptr<GammaSurface> gammaSurface(const size_t&) const;

    };
}
#endif
