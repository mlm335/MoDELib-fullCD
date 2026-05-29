/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SlipSystemTab_H_
#define model_SlipSystemTab_H_

#include <map>
#include <memory>

#include <QWidget>
#include <QGridLayout>
#include <QPushButton>


#include <DefectiveCrystal.h>
#include <SlipSystem.h>

namespace model
{

class SlipSystemButton : public QPushButton
 {
     Q_OBJECT // Essential for signals/slots

 public:
     
     const int gID;
     const int sID;
     
     SlipSystemButton(const int& gID_in, const int& sID_in,const QString &text, QWidget *parent);
     virtual ~SlipSystemButton();
     Eigen::Matrix<int,3,1> rgb() const;
     QColor qColor() const;

 };

    struct SlipSystemTab : public QWidget
    {
        
        Q_OBJECT
        
        
    private slots:
        
        void changeColor();
        void setBurgersColor();
        void setNormalColor();
        void setFamilyColor();
        
//        void modify();
//        void reload();
//        void updatePlots();

        
    public:

        typedef std::map<const SlipSystem*,SlipSystemButton*> SlipSystemColorMapType;
        
        const int numCols;
        SlipSystemColorMapType slipSystemColorMap;
        QGridLayout* mainLayout;

        

        SlipSystemTab(const  DefectiveCrystal<3>&);
        
//        void updateConfiguration(const size_t& frameID);

        
    };
    
} // namespace model
#endif







