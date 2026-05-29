/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InterpolantBase_h_
#define model_InterpolantBase_h_

#include <map>
#include <Eigen/Dense>

namespace model
{

    struct ExtrapolationMethod
    {
        int type;
        double period;
        
        ExtrapolationMethod();
        ExtrapolationMethod(const int&,const double&);

    };
    
    class InterpolantBase : public std::map<double,double>
    {

    public:
        typedef std::map<double,double> MapType;
        typedef Eigen::Matrix<double,Eigen::Dynamic,2> MatrixType;
        ExtrapolationMethod extrapolation; // 0=default, 1=constant

    private:
        static MapType matrix2Map(const MatrixType& data);
        static MatrixType map2Matrix(const MapType& data);
        
        virtual double _f(const double& x) const = 0;


    public:
                
        InterpolantBase(const MapType& data,const ExtrapolationMethod& extrap);
        InterpolantBase(const MatrixType& data,const ExtrapolationMethod& extrap);
        virtual ~InterpolantBase();
        double xMin() const;
        double xMax() const;
        double operator()(const double& x) const; // pure virtual
        double f(const double& x) const;
        double atPeriodic(const double& x,const double& xL,const double& xH) const;
        MapType& map();
        const MapType& map() const;

    };
    
} // close namespace model
#endif
