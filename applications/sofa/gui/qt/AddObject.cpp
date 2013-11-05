/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#include "AddObject.h"
#include "RealGUI.h"
#include "FileManagement.h"

#include <iostream>
#include <sstream>

#ifdef SOFA_QT4
#include <Q3FileDialog>
#include <QLineEdit>
#include <QLabel>
#include <QRadioButton>
#include <QPushButton>
#include <Q3ButtonGroup>
#include <QGridLayout>
#include <qevent.h>
#else
#include <qfiledialog.h>
#include <qlineedit.h>
#include <qlabel.h>
#include <qradiobutton.h>
#include <qpushbutton.h>
#include <qbuttongroup.h>
#include <qlayout.h>
#endif

namespace sofa
{

namespace gui
{

namespace qt
{


#ifndef SOFA_QT4
typedef QFileDialog  Q3FileDialog;
typedef QButtonGroup Q3ButtonGroup;
#endif


  AddObject::AddObject( std::vector< std::string > *list_object_, QWidget* parent, bool , Qt::WFlags ): list_object(list_object_)
{
    setupUi(this);
    //At the creation of the dialog window, we enable the custom object
    custom->setChecked(true);


    //Creation of the list of radio button corresponding to the preset objects: they are specified in the sofa/scenes/object.txt file
    if (list_object != NULL)
    {
        QRadioButton *button;
        std::string current_name;

        for (int i=0; i<(int)list_object->size(); i++)
        {
            std::ostringstream ofilename;
            current_name = (*list_object)[i];
            std::string::size_type pos=current_name.rfind('/');

            if (pos != std::string::npos)
            {
                current_name = current_name.substr(pos+1, current_name.size()-pos-5);
            }
            button = new QRadioButton( buttonGroup, QString(current_name.c_str()) );
            button->setText(current_name.c_str());

#ifdef SOFA_QT4
            gridLayout1->addWidget( button, i+1, 0 );
#else
            buttonGroupLayout->addWidget( button, i+1, 0 );
#endif
        }
    }
    positionX->setText("0");
    positionY->setText("0");
    positionZ->setText("0");

    rotationX->setText("0");
    rotationY->setText("0");
    rotationZ->setText("0");

    scaleValue->setText("1");
    //Option still experimental : disabled
    scaleValue->hide();
    scaleText->hide();

    openFilePath->setText(NULL);

    //Make the connection between this widget and the parent
    connect( this, SIGNAL(loadObject(std::string, double, double, double, double, double, double,double)), parent, SLOT(loadObject(std::string, double, double, double,double, double, double, double)));
    //For tje Modifications of the state of the radio buttons
    connect( buttonGroup, SIGNAL( clicked(int) ), this, SLOT (buttonUpdate(int)));
}

//**************************************************************************************
//When the Ok Button is clicked, this method is called: we just have to emit a signal to the parent, with the information on the object
void AddObject::accept()
{
    std::string position[3];
    std::string rotation[3];
    std::string scale;
#ifdef SOFA_QT4
    std::string object_fileName(openFilePath->text().toStdString());
    position[0] = positionX->text().toStdString();
    position[1] = positionY->text().toStdString();
    position[2] = positionZ->text().toStdString();

    rotation[0] = rotationX->text().toStdString();
    rotation[1] = rotationY->text().toStdString();
    rotation[2] = rotationZ->text().toStdString();

    scale       = scaleValue->text().toStdString();
#else
    std::string object_fileName(openFilePath->text().latin1());
    position[0] = positionX->text().latin1();
    position[1] = positionY->text().latin1();
    position[2] = positionZ->text().latin1();

    rotation[0] = rotationX->text().latin1();
    rotation[1] = rotationY->text().latin1();
    rotation[2] = rotationZ->text().latin1();

    scale       = scaleValue->text().latin1();
#endif
    emit( loadObject(object_fileName, atof(position[0].c_str()),atof(position[1].c_str()),atof(position[2].c_str()),
            atof(rotation[0].c_str()),atof(rotation[1].c_str()),atof(rotation[2].c_str()),
            atof(scale.c_str())));
    setPath(object_fileName);
    QDialog::accept();
}

//**************************************************************************************
//Set the default file
void AddObject::setPath(const std::string path)
{
    fileName = path;
    openFilePath->setText(QString(fileName.c_str()));
}

//**************************************************************************************
//Open a file Dialog and set the path of the selected path in the text field.
void AddObject::fileOpen()
{
    QString s  = getOpenFileName(this, QString(fileName.c_str()), "Scenes (*.xml *.scn);;All (*)", "open file dialog",  "Choose a file to open" );

    if (s.isNull() ) return;
#ifdef SOFA_QT4
    std::string object_fileName(s.toStdString());
#else
    std::string object_fileName(s.latin1());
#endif

    openFilePath->setText(QString(object_fileName.c_str()));
}

//**************************************************************************************
//The state of the radio buttons has been modified
//we update the content of the dialog window
void AddObject::buttonUpdate(int Id)
{
    //Id = 0 : custom radio button clicked: we need to show the selector of the file
    if (Id == 0)
    {
        openFilePath->setText(fileName.c_str());
        openFileText->show();
        openFilePath->show();
        openFileButton->show();
    }
    else
    {
        openFilePath->setText((*list_object)[Id-1].c_str());
        openFileText->hide();
        openFilePath->hide();
        openFileButton->hide();
    }
}
} // namespace qt

} // namespace gui

} // namespace sofa

