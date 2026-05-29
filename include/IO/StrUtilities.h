/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StrUtilities_H_
#define model_StrUtilities_H_

//#include <filesystem>
#include <string>
#include <algorithm> // Required for std::transform
#include <cctype>    // Required for std::tolower

namespace model
{

    struct StrUtilities
    {
      
        static std::string removeSpaces(std::string key)
        {
            key.erase(std::remove_if(key.begin(), key.end(), [](unsigned char x) { return std::isspace(x); }), key.end());
            return key;
        }
        
        static std::string lowercase(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            return text;
        }
        
    };

}
#endif
