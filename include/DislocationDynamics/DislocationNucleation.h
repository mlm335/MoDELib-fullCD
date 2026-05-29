/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationNucleation_H
#define model_DislocationNucleation_H

#include <deque>
#include <map>
#include <tuple>

#include <Eigen/Dense>

//#include <TypeTraits.h>
//#include <LatticeModule.h>
//#include <BCClattice.h>
//#include <FCClattice.h>
//#include <Grain.h>
//#include <TerminalColors.h>
#include <DislocationNetwork.h>

namespace model
{

    template <int dim>
    struct DislocationNucleation
    {
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        
        DislocationNetwork<dim>& DN;
        
        DislocationNucleation(DislocationNetwork<dim>& DN_in);
        
        void bulkNucleate(const int& bulkNucleationModel);
        void surfaceNucleate(const int& surfaceNucleationModel);
    };

}
#endif
