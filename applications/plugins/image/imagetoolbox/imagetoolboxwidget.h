/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/


#ifndef IMAGETOOLBOXWIDGET_H
#define IMAGETOOLBOXWIDGET_H


#include "imagetoolboxcentralwidget.h"

#include "imagetoolboxbasicactionwidget.h"
#include "imagetoolboxlabelactionwidget.h"

namespace sofa
{
namespace gui
{
namespace qt
{

using namespace helper;
using namespace cimg_library;
using defaulttype::Vec;


class imagetoolbox_data_widget_container: public QObject
{
    Q_OBJECT

public slots:
    virtual void handleSliderPolicies()=0; // needed for synchronization of slider visiblity
    //virtual void draw()=0;

signals:
    //void updateImage();

};

template<class T>
class Timagetoolbox_data_widget_container : public imagetoolbox_data_widget_container
{


    typedef T ImageToolBoxDataType;
    typedef TImageToolBoxCentralWidget<T> CentralWidget;
    typedef ImageToolBoxBasicActionWidget BasicToolBar;
    typedef ImageToolBoxLabelActionWidget LabelToolBar;
    
    typedef typename ImageToolBoxDataType::ImagePlaneType ImagePlaneType;
//    typedef TImagePlaneGraphWidget<ImagePlaneType> Graph;
    
public:
    
    typedef QVBoxLayout Layout;

    Layout* container_layout;
    
    QMainWindow* main;
    CentralWidget* central;
    BasicToolBar* basicTB;
    LabelToolBar* labelTB;
    

    Timagetoolbox_data_widget_container() : /*setting(NULL), graph(NULL), options(NULL),*/ container_layout(NULL),main(NULL),basicTB(NULL),labelTB(NULL) {}

    bool createLayout( DataWidget* parent )
    {
        if( parent->layout() != NULL || container_layout != NULL ) return false;
        container_layout = new Layout(parent);
        return true;
    }

    bool createLayout( QLayout* layout)
    {
        if ( container_layout != NULL ) return false;
        container_layout = new Layout(layout);
        return true;
    }

    bool createWidgets(DataWidget* /*parent*/, const ImageToolBoxDataType& d, bool /*readOnly*/)
    {
        //std::cout << "imagetoolbox_data_widget_container::createWidgets" << std::endl;
        main = new QMainWindow();
        
        // std::cout << "a" << main << std::endl;
        central = new CentralWidget();
        
         //std::cout << "b" << central << std::endl;
        central->createLayout();
        central->createWidgets(d.getPlane());
        central->insertWidgets();
        
        basicTB = new BasicToolBar();
        basicTB->connectCentralW(central);
        
        labelTB = new LabelToolBar();
        labelTB->connectCentralW(central);
        labelTB->setGraphScene(central->graphXY->graphscene(),central->graphXZ->graphscene(),central->graphZY->graphscene());
        labelTB->setLabels(d.getLabels());
        //connect(labelTB,SIGNAL(updateImage()),this,SIGNAL(updateImage()));
        //connect(this,SIGNAL(updateImage()),this,SLOT(draw()));;
        //connect(this,SIGNAL(updateImage()),central,SLOT(update()));
        //connect(this,SIGNAL(updateImage()),this,SLOT(testsignal()));

        //std::cout << "c" << std::endl;
        
        main->setCentralWidget(central);//central);
        
        main->addToolBar(Qt::TopToolBarArea,basicTB);

        QDockWidget *dw = new QDockWidget("Label Tools");
        dw->setWidget(labelTB);
        main->addDockWidget(Qt::RightDockWidgetArea,dw);
        
    /*
        setting = new Setting(parent);
        setting->readFromData(d);


        graph = new Graph(setting,parent);

        if(setting->getChannelMax())
        {
            options = new Options(setting,graph,parent);
            options->setRange(0,setting->getChannelMax());
        }
        else options=NULL;*/

        //std::cout << "~imagetoolbox_data_widget_container::createWidgets" << std::endl;
        return true;
    }

    void setReadOnly(bool /*readOnly*/) { }
    void readFromData(const ImageToolBoxDataType& d0)
    { 
        if(central) central->readFromData(d0.getPlane());
        
        //const helper::vector< sofa::component::engine::LabelImageToolBox*> &dd = d0.;
        
        if(labelTB)
        {
            //labelTB->setLabels(d0.getLabels());
        }
        
         /*if(graphXY) graphXY->readFromData(d0);
         if(graphXZ) graphXZ->readFromData(d0);
         if(graphZY) graphZY->readFromData(d0);*/
    }
    
    void writeToData(ImageToolBoxDataType& d)
    {
        if(central) central->writeToData(d.plane());
        
       /* if(graphXY) graphXY->writeToData(d);
        if(graphXZ) graphXZ->writeToData(d);
        if(graphZY) graphZY->writeToData(d);*/
    }
    
    void insertWidgets()
    {
        assert(container_layout);
        //if(label) container_layout->add(label);
        
        if(main)container_layout->add(main);
        
       // if(setting) container_layout->add(setting->getWidget());
       // if(options) container_layout->add(options);
    }
    
        void handleSliderPolicies()
        {
            central->handleSliderPolicies();
        }

        /*virtual void draw()
        {
            std::cout << "draw" <<std::endl;
            this->central->graphXY->draw();
            this->central->graphZY->draw();
            this->central->graphXZ->draw();
        }*/
};


template<class T>
class SOFA_IMAGE_API ImageToolBoxWidget : public SimpleDataWidget<T, Timagetoolbox_data_widget_container< T > >
{
    typedef SimpleDataWidget<T, Timagetoolbox_data_widget_container< T > > Inherit;
    typedef sofa::core::objectmodel::Data<T> MyData;
public:
    ImageToolBoxWidget(QWidget* parent,const char* name, MyData* d) : Inherit(parent,name,d) {}
    
    virtual unsigned int sizeWidget() {return 3;}
    virtual unsigned int numColumnWidget() {return 1;}
    virtual bool createWidgets()
    {
    
    //std::cout << "ImageToolBoxWidget::createWidgets" << std::endl;
        bool b = Inherit::createWidgets();
        
    //std::cout << "~ImageToolBoxWidget::createWidgets" << std::endl;
    
        imagetoolbox_data_widget_container* s = dynamic_cast<imagetoolbox_data_widget_container *>(&this->container);
        //this->connect(s,SIGNAL(updateImage()),this,SLOT(setWidgetDirty()));

        //this->connect(s,SIGNAL(updateImage()),this,SLOT(forceUpdateWidgetValue()));
        //this->connect(s,SIGNAL(updateImage()),this,SLOT(updateDataValue()));




        //this->connect(s,SIGNAL(clampModified()), this, SLOT(setWidgetDirty()));
        return b;
    }
};

}

}

}


#endif // IMAGETOOLBOXWIDGET_H
