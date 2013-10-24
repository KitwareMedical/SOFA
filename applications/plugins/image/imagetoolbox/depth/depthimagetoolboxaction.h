#ifndef DEPTHIMAGETOOLBOXACTION_H
#define DEPTHIMAGETOOLBOXACTION_H

#include <QtGui>


#include "../labelimagetoolboxaction.h"
#include "../labelgrid/labelgridimagetoolbox.h"
//#include "depthimagetoolbox.h"

#include "initImage.h"

namespace sofa
{

namespace component
{

namespace engine
{
class DepthImageToolBox;
}}}

namespace sofa
{
namespace gui
{
namespace qt
{

class SOFA_IMAGE_API DepthRowImageToolBoxAction:public QObject
{
    Q_OBJECT

    QComboBox *gridchoice1, *gridchoice2, *basechoice;
    QLineEdit *offset1, *offset2;
    QTableWidgetItem *name;
    QWidget *w1,*w2;

    int index;
    static int numindex;

public:
    DepthRowImageToolBoxAction():QObject(),index(0)
    {
        gridchoice1 = new QComboBox();
        offset1 = new QLineEdit();
        QFormLayout *layout1 = new QFormLayout();
        layout1->addRow("grid",gridchoice1);
        layout1->addRow("offset",offset1);
        w1=new QWidget(); w1->setLayout(layout1);

        gridchoice2 = new QComboBox();
        offset2 = new QLineEdit();
        QFormLayout *layout2 = new QFormLayout();
        layout2->addRow("grid",gridchoice2);
        layout2->addRow("offset",offset2);
        w2=new QWidget(); w2->setLayout(layout2);

        basechoice = new QComboBox();

        name = new QTableWidgetItem(QString("layer") + QString::number(numindex));

        connectChange();
    }

    ~DepthRowImageToolBoxAction()
    {
        delete w1;
        delete w2;
        delete name;
    }

    void setParameters(helper::vector< sofa::component::engine::LabelGridImageToolBoxNoTemplated* > &list)
    {


        for(unsigned int i=0;i<list.size();i++)
        {
            QString str(list[i]->getName().c_str());
            gridchoice1->addItem(str);
            gridchoice2->addItem(str);
            basechoice->addItem(str);
        }
        offset1->setText("0");
        offset2->setText("0");
    }
    
    void setValue(QString _name, int _layerId1, QString _offset1, int _layerId2, QString _offset2,int _base)
    {
        disconnectChange();
        name->setText(_name);
        gridchoice1->setCurrentIndex(_layerId1);
        gridchoice2->setCurrentIndex(_layerId2);
        basechoice->setCurrentIndex(_base);
        offset1->setText(_offset1);
        offset2->setText(_offset2);
        connectChange();
        emit valueChanged(index,name->text(),gridchoice1->currentIndex(),offset1->text(),gridchoice2->currentIndex(),offset2->text(),basechoice->currentIndex());
    }


    void toTableWidgetRow(QTableWidget *listLayers)
    {
        index = listLayers->rowCount();
        std::cout << "index"<<index<<std::endl;


        listLayers->insertRow(index);

        w1->adjustSize();
        w2->adjustSize();
        basechoice->adjustSize();

        int size = w1->minimumHeight();if(size<w2->minimumHeight())size=w2->minimumHeight();


        listLayers->setItem(index,0,name);
        listLayers->setCellWidget(index,1,basechoice);
        listLayers->setCellWidget(index,2,w1);
        listLayers->setCellWidget(index,3,w2);

        listLayers->setRowHeight(index,size);


        int size1 = basechoice->minimumWidth();
        if(size1>listLayers->columnWidth(1))listLayers->setColumnWidth(1,size1);
        int size2 = w1->minimumWidth();
        if(size2>listLayers->columnWidth(2))listLayers->setColumnWidth(2,size2);
        int size3 = w2->minimumWidth();
        if(size3>listLayers->columnWidth(3))listLayers->setColumnWidth(3,size3);


    }
    
private:
    inline void connectChange()
    {
        this->connect(gridchoice1,SIGNAL(currentIndexChanged(int)),this,SLOT(change()));
        this->connect(gridchoice2,SIGNAL(currentIndexChanged(int)),this,SLOT(change()));
        this->connect(basechoice,SIGNAL(currentIndexChanged(int)),this,SLOT(change()));
        this->connect(offset1,SIGNAL(textEdited(QString)),this,SLOT(change()));
        this->connect(offset2,SIGNAL(textEdited(QString)),this,SLOT(change()));
    }
    
    inline void disconnectChange()
    {
        this->disconnect(gridchoice1,SIGNAL(currentIndexChanged(int)),this,SLOT(change()));
        this->disconnect(gridchoice2,SIGNAL(currentIndexChanged(int)),this,SLOT(change()));
        this->disconnect(basechoice,SIGNAL(currentIndexChanged(int)),this,SLOT(change()));
        this->disconnect(offset1,SIGNAL(textEdited(QString)),this,SLOT(change()));
        this->disconnect(offset2,SIGNAL(textEdited(QString)),this,SLOT(change()));
    }

public slots:
    void change()
    {
        emit valueChanged(index,name->text(),gridchoice1->currentIndex(),offset1->text(),gridchoice2->currentIndex(),offset2->text(),basechoice->currentIndex());
    }

signals:
    void valueChanged(int id, QString name, int layerId1, QString offset1, int layerId2, QString offset2,int);
};


class SOFA_IMAGE_API DepthImageToolBoxAction : public LabelImageToolBoxAction
{
Q_OBJECT



    //QGraphicsLineItem *lineH[3], *lineV[3];
    QGraphicsPathItem *path[3];

public:



    typedef sofa::component::engine::DepthImageToolBox DepthImageToolBox;
    typedef sofa::component::engine::LabelGridImageToolBoxNoTemplated LabelGridImageToolBoxNoTemplated;



public:
    DepthImageToolBoxAction(sofa::component::engine::LabelImageToolBox* lba,QObject *parent);
    ~DepthImageToolBoxAction();
    
    sofa::component::engine::DepthImageToolBox* DITB();

private:
    void createMainCommands();
    void createListLayers();
    void textToOffset(QString text, double &outValue, int &type);
    void offsetToText(QString &out, double outValue, int type);


    void lineTo(const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition);
    void moveTo(const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition);


public slots:
    virtual void addOnGraphs();
    virtual void updateGraphs();
    virtual void updateColor();



private slots:
    void selectionPointEvent(int mouseevent, const unsigned int axis,const sofa::defaulttype::Vec3d& imageposition,const sofa::defaulttype::Vec3d& position3D,const QString& value);
    
    void executeButtonClick();
    void saveButtonClick();
    void saveSceneButtonClick();
    void loaduttonClick();

    void createNewRow();
    void createNewRow(int layer);
    void changeRow(int, QString, int, QString, int, QString, int base);
    
private:
    QPushButton* select;
    QTableWidget *listLayers;
    helper::vector< DepthRowImageToolBoxAction* > listRows;
    
};

}
}
}

#endif // DepthImageToolBoxACTION_H



