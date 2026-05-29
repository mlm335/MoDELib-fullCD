/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SlipSystemTab_cpp_
#define model_SlipSystemTab_cpp_

#include <QColorDialog>
#include <QPalette>
#include <QColor>

#include <vtkLookupTable.h>

#include <SlipSystemTab.h>
#include <Vector2Color.h>



namespace model
{

     SlipSystemButton::SlipSystemButton(const int& gID_in, const int& sID_in,const QString &text, QWidget *parent = nullptr):
     /* init */ QPushButton(text,parent)
     /* init */,gID(gID_in)
     /* init */,sID(sID_in)
     {
         
     }

    SlipSystemButton::~SlipSystemButton()
    {
        
    }

Eigen::Matrix<int,3,1> SlipSystemButton::rgb() const
{
    QColor backgroundColor(qColor());
    return (Eigen::Matrix<int,3,1>()<<backgroundColor.red(), backgroundColor.green(), backgroundColor.blue()).finished();
}

QColor SlipSystemButton::qColor() const
{
    return this->palette().color(QPalette::Window);
}




SlipSystemTab::SlipSystemTab(const  DefectiveCrystal<3>& defectiveCrystal) :
/* init */ numCols(4)
/* init */,mainLayout(new QGridLayout(this))
        {
            
            QPushButton* burgersButton(new QPushButton("slipSystem Burgers"));
            connect( burgersButton, SIGNAL(clicked()), this, SLOT(setBurgersColor()) );
            mainLayout->addWidget(burgersButton,0,0,1,numCols);

            QPushButton* normalButton(new QPushButton("slipSystem normal"));
            connect( normalButton, SIGNAL(clicked()), this, SLOT(setNormalColor()) );
            mainLayout->addWidget(normalButton,1,0,1,numCols);

            QPushButton* familyButton(new QPushButton("slipSystem family"));
            connect( familyButton, SIGNAL(clicked()), this, SLOT(setFamilyColor()) );
            mainLayout->addWidget(familyButton,2,0,1,numCols);

            
            size_t ssCount(0);
            for(const auto& grain : defectiveCrystal.ddBase.poly.grains)
            {
                for(const auto& slipSystem : grain.second->slipSystems())
                {
//                    const std::pair<int,int> key(std::make_pair(grain.first,slipSystem.second->sID));
                    
                    const int i = ssCount / numCols+3;
                    const int j = ssCount % numCols;
                    
                    SlipSystemButton* ssb(new SlipSystemButton(grain.first,slipSystem.second->sID,QString::fromStdString("g"+std::to_string(grain.first)+",s"+std::to_string(slipSystem.second->sID))));
                    connect( ssb, SIGNAL(clicked()), this, SLOT(changeColor()) );
                    slipSystemColorMap.emplace(slipSystem.second.get(),ssb);
                    mainLayout->addWidget(ssb,i,j,1,1);
                    ssCount++;
                }
            }
            this->setLayout(mainLayout);
        }
        
void SlipSystemTab::changeColor()
{
    SlipSystemButton* clickedButton(qobject_cast<SlipSystemButton*>(sender()));
    if(clickedButton)
    {
        QColor newColor = QColorDialog::getColor(clickedButton->qColor());
        clickedButton->setStyleSheet(QString::fromStdString("background-color: rgb("+std::to_string(newColor.red())+","
                                                            /*                    */+std::to_string(newColor.green())+","
                                                            /*                    */+std::to_string(newColor.blue())+");"));
        
    }

//    std::cout<< key.first<<" "<< key.second<<std::endl;
//    if ( newColor != color )
//    {
//        setColor( newColor );
//    }
}

void SlipSystemTab::setBurgersColor()
{
    for(auto& pair : slipSystemColorMap)
    {
        //const auto clr(Vector2Color::v2c(pair.first->unitSlip));
        const auto clr(Vector2Color::v2c(pair.first->s.dir.base().cast<double>()));
        pair.second->setStyleSheet(QString::fromStdString("background-color: rgb("+std::to_string(clr(0))+","
                                                            /*                    */+std::to_string(clr(1))+","
                                                            /*                    */+std::to_string(clr(2))+");"));

    }
}

void SlipSystemTab::setNormalColor()
{
    for(auto& pair : slipSystemColorMap)
    {
//        const auto clr(Vector2Color::v2c(pair.first->unitNormal));
        const auto clr(Vector2Color::v2c(pair.first->n.base().cast<double>()));
        pair.second->setStyleSheet(QString::fromStdString("background-color: rgb("+std::to_string(clr(0))+","
                                                            /*                    */+std::to_string(clr(1))+","
                                                            /*                    */+std::to_string(clr(2))+");"));

    }
}

void SlipSystemTab::setFamilyColor()
{
    
    std::map<std::pair<float,float>,std::set<const SlipSystem*>> ssFamilyMap;
    
    for(auto& pair : slipSystemColorMap)
    {
        const std::pair<float,float> key(pair.first->s.cartesian().squaredNorm(),pair.first->n.cartesian().squaredNorm());
        ssFamilyMap[key].insert(pair.first);
    }
    
    vtkSmartPointer<vtkLookupTable> familyLut(vtkSmartPointer<vtkLookupTable>::New());
//    familyLut->SetHueRange(0.66667, 0.0);
    familyLut->SetTableRange(0, ssFamilyMap.size());
    familyLut->Build();

    int fID(0);
    for(const auto& pair : ssFamilyMap)
    {
        for(const auto& ss : pair.second)
        {
            const auto ssIter(slipSystemColorMap.find(ss));
            if(ssIter!=slipSystemColorMap.end())
            {
                double dclr[3];
                familyLut->GetColor(fID, dclr);
                ssIter->second->setStyleSheet(QString::fromStdString("background-color: rgb("+std::to_string(int(dclr[0]*255))+","
                                                                    /*                    */+std::to_string(int(dclr[1]*255))+","
                                                                    /*                    */+std::to_string(int(dclr[2]*255))+");"));
            }
        }
        fID++;
    }
}


} // namespace model
#endif







