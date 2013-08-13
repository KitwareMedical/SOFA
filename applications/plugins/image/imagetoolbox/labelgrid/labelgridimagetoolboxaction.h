#ifndef LABELGRIDIMAGETOOLBOXACTION_H
#define LABELGRIDIMAGETOOLBOXACTION_H

#include <QtGui>


#include "../labelimagetoolboxaction.h"
//#include "labelpointimagetoolbox.h"

#include "initImage.h"

namespace sofa
{

namespace component
{

namespace engine
{
class LabelGridImageToolBoxNoTemplated;
}}}

namespace sofa
{
namespace gui
{
namespace qt
{

class SOFA_IMAGE_API LabelGridImageToolBoxAction : public LabelImageToolBoxAction
{
Q_OBJECT

    //QGraphicsLineItem *lineH[3], *lineV[3];
    QGraphicsPathItem *path[3];

public:
    LabelGridImageToolBoxAction(sofa::component::engine::LabelImageToolBox* lba,QObject *parent);
    ~LabelGridImageToolBoxAction();
    
    sofa::component::engine::LabelGridImageToolBoxNoTemplated* LGITB();

private:
    void createMainCommands();
    //void createBoxCommands();
    void createGridCommands();
    //void createNormalCommands();

    void lineTo(const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition);
    void moveTo(const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition);


public slots:
    virtual void addOnGraphs();
    virtual void updateGraphs();
    virtual void updateColor();

private slots:
    void selectionPointEvent(int mouseevent, const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition,const sofa::defaulttype::Vec3d& position3D,const QString& value);
    
    void executeButtonClick();
    
private:
    QPushButton* select;
    
};

}
}
}

#endif // LabelGridIMAGETOOLBOXACTION_H



