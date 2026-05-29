/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearInterpolant_H_
#define model_LinearInterpolant_H_

#include <InterpolantBase.h>

namespace model{
    
    class LinearInterpolant : public InterpolantBase
    {
        
        static double linearSegment(const InterpolantBase::MapType::const_iterator& iter1,const InterpolantBase::MapType::const_iterator& iter2,const double& x);
        double _f(const double& x) const override;

    public:
        
        LinearInterpolant(const InterpolantBase::MatrixType& data,const ExtrapolationMethod& extrap=ExtrapolationMethod());
        LinearInterpolant(const InterpolantBase::MapType& data,const ExtrapolationMethod& extrap=ExtrapolationMethod());

    };
    
} // close namespace model
#endif
