#ifndef LABELPOINTIMAGETOOLBOXACTION_H
#define LABELPOINTIMAGETOOLBOXACTION_H

#include <QAction>
#include <QGraphicsLineItem>


#include "../labelimagetoolboxaction.h"
//#include "labelpointimagetoolbox.h"

#include "initImage.h"

namespace sofa
{

namespace component
{

namespace engine
{
class LabelPointImageToolBox;
}}}

namespace sofa
{
namespace gui
{
namespace qt
{

class SOFA_IMAGE_API LabelPointImageToolBoxAction : public LabelImageToolBoxAction
{
Q_OBJECT

    QGraphicsLineItem *lineH[3], *lineV[3];

public:
    LabelPointImageToolBoxAction(sofa::component::engine::LabelImageToolBox* lba,QObject *parent);
    ~LabelPointImageToolBoxAction();
    
    sofa::component::engine::LabelPointImageToolBox* LPITB();

private:

public slots:
    virtual void addOnGraphs();
    virtual void updateGraphs();
    virtual void updateColor();
    
private slots:
    void selectionPointButtonClick(bool);
    void selectionPointEvent(int mouseevent, const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition,const sofa::defaulttype::Vec3d& position3D,const QString& value);
    void sectionButtonClick();
    
    
private:
    QAction* select;
    
};

}
}
}

#endif // LABELPOINTIMAGETOOLBOXACTION_H
