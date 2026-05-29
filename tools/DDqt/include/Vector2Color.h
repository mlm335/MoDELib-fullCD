/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Vector2Color_H_
#define model_Vector2Color_H_


#include <cfloat>
#include <Eigen/Dense>

namespace model
{
    struct Vector2Color
    {
        
        static Eigen::Matrix<int,3,1> v2c(Eigen::Matrix<double,3,1> clrVector)
        {
            float clrTol=100.0*FLT_EPSILON;
            if(clrVector(0)<-clrTol)
            {// first component not zero but begative, flip color
                clrVector*=-1.0;
            }
            else if(fabs(clrVector(0))<=clrTol)
            {// first component is zero, use second component
                if(clrVector(1)<-clrTol)
                {// second component not zero but begative, flip color
                    clrVector*=-1.0;
                }
                else if(fabs(clrVector(1))<=clrTol)
                {// second component is zero, use third component
                    if(clrVector(2)<-clrTol)
                    {
                        clrVector*=-1.0;
                    }
                }
            }
            
            clrVector = (clrVector + Eigen::Matrix<double,3,1>::Ones(3) * clrVector.norm()).eval();
            clrVector.normalize();
            return (clrVector*255).cast<int>();
        }
        
    };
}
#endif







