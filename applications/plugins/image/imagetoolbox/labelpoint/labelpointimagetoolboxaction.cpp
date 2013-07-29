#include <QString>

#include "labelpointimagetoolboxaction.h"

#include "labelpointimagetoolbox.h"

namespace sofa
{
namespace gui
{
namespace qt
{

LabelPointImageToolBoxAction::LabelPointImageToolBoxAction(sofa::component::engine::LabelImageToolBox* lba,QObject *parent):
    LabelImageToolBoxAction(lba,parent)
{
    //button selection point
    select = new QAction(this);
    this->l_actions.append(select);
    select->setText("Select Point");
    select->setCheckable(true);
    connect(select,SIGNAL(toggled(bool)),this,SLOT(selectionPointButtonClick(bool)));
    
    QAction* section = new QAction(this);
    this->l_actions.append(section);
    section->setText("Section");
    connect(section,SIGNAL(triggered()),this,SLOT(sectionButtonClick()));
    
}

LabelPointImageToolBoxAction::~LabelPointImageToolBoxAction()
{
    //delete select;
}


sofa::component::engine::LabelPointImageToolBox* LabelPointImageToolBoxAction::LPITB()
{
    return dynamic_cast<sofa::component::engine::LabelPointImageToolBox*>(this->p_label);
}

void LabelPointImageToolBoxAction::selectionPointEvent(int /*mouseevent*/, const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition,const sofa::defaulttype::Vec3d& position3D,const QString& value)
{
    
    select->setChecked(false);
    disconnect(this,SIGNAL(clickImage(int,uint,sofa::defaulttype::Vec3d,sofa::defaulttype::Vec3d,QString)),this,SLOT(selectionPointEvent(int,uint,sofa::defaulttype::Vec3d,sofa::defaulttype::Vec3d,QString)));
    
    sofa::component::engine::LabelPointImageToolBox* lp = LPITB();
    
    lp->d_ip.setValue(imageposition);
    lp->d_p.setValue(position3D);
    lp->d_axis.setValue(axis);
    lp->d_value.setValue(value.toStdString());
    
    updateGraphs();
}

void LabelPointImageToolBoxAction::selectionPointButtonClick(bool b)
{
    
    if(b)
    {
        //select->setChecked(true);
        connect(this,SIGNAL(clickImage(int,uint,sofa::defaulttype::Vec3d,sofa::defaulttype::Vec3d,QString)),this,SLOT(selectionPointEvent(int,uint,sofa::defaulttype::Vec3d,sofa::defaulttype::Vec3d,QString)));
    }
    else
    {
        //select->setChecked(false);
        disconnect(this,SIGNAL(clickImage(int,uint,sofa::defaulttype::Vec3d,sofa::defaulttype::Vec3d,QString)),this,SLOT(selectionPointEvent(int,uint,sofa::defaulttype::Vec3d,sofa::defaulttype::Vec3d,QString)));
    }
    
}

void LabelPointImageToolBoxAction::addOnGraphs()
{

//    std::cout << "addOnGraph"<<std::endl;

    lineH[0] = GraphXY->addLine(0,0,0,0);
    lineH[1] = GraphXZ->addLine(0,0,0,0);
    lineH[2] = GraphZY->addLine(0,0,0,0);
    
    lineV[0] = GraphXY->addLine(0,0,0,0);
    lineV[1] = GraphXZ->addLine(0,0,0,0);
    lineV[2] = GraphZY->addLine(0,0,0,0);
    
    for(int i=0;i<3;i++)
    {
        lineH[i]->setVisible(false);
        lineV[i]->setVisible(false);
    }
    
    updateColor();
}

void LabelPointImageToolBoxAction::updateGraphs()
{
    sofa::defaulttype::Vec3d pos = LPITB()->d_ip.getValue();
    
    QRectF boundaryXY = GraphXY->itemsBoundingRect();
    
//    std::cout << "updateOnGraphs"<<std::endl;
    
    lineH[0]->setVisible(true);
    lineH[0]->setLine(pos.x()-4,pos.y(),pos.x()+4,pos.y());
    
    lineV[0]->setVisible(true);
    lineV[0]->setLine(pos.x(),pos.y()-4,pos.x(),pos.y()+4);
    
    
    lineH[1]->setVisible(true);
    lineH[1]->setLine(pos.x()-4,pos.z(),pos.x()+4,pos.z());
    
    lineV[1]->setVisible(true);
    lineV[1]->setLine(pos.x(),pos.z()-4,pos.x(),pos.z()+4);
    
    
    lineH[2]->setVisible(true);
    lineH[2]->setLine(pos.z()-4,pos.y(),pos.z()+4,pos.y());
    
    lineV[2]->setVisible(true);
    lineV[2]->setLine(pos.z(),pos.y()-4,pos.z(),pos.y()+4);
    
}

void LabelPointImageToolBoxAction::updateColor()
{
    for(int i=0;i<3;i++)
    {
        lineH[i]->setPen(QPen(this->color()));
        lineV[i]->setPen(QPen(this->color()));
    }
}

void LabelPointImageToolBoxAction::sectionButtonClick()
{
   // std::cout << "LabelPointImageToolBoxAction::sectionButtonClick()"<<std::endl;
    sofa::defaulttype::Vec3d pos = LPITB()->d_ip.getValue();
    
    sofa::defaulttype::Vec3i pos2(round(pos.x()),round(pos.y()),round(pos.z()));

    emit sectionChangeGui(pos2);
}












SOFA_DECL_CLASS(LabelPointImageToolBoxAction)



}
}
}
