/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LagrangeInterpolant_H_
#define model_LagrangeInterpolant_H_

#include <InterpolantBase.h>

namespace model
{

    class LagrangeInterpolant : public InterpolantBase
    {
        double p(const InterpolantBase::const_iterator&,const double& x) const;
        double lagrangePoly(const double& x) const;
        double _f(const double& x) const override;

    public:
        
        LagrangeInterpolant(const InterpolantBase::MatrixType& data,const ExtrapolationMethod& extrap=ExtrapolationMethod());
        LagrangeInterpolant(const InterpolantBase::MapType& data,const ExtrapolationMethod& extrap=ExtrapolationMethod());

	};
    
}
#endif
