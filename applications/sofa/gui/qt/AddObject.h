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
#ifndef ADDOBJECT_H
#define ADDOBJECT_H


#include "ui_DialogAddObject.h"
#include <vector>

namespace sofa
{

namespace gui
{

namespace qt
{



class AddObject : public QDialog, public Ui_DialogAddObject
{
    Q_OBJECT
public:

    AddObject( std::vector<std::string> *list_object_, QWidget* parent, bool  modal= FALSE, Qt::WFlags f= 0 );


    void setPath(const std::string path);

public slots:
    void fileOpen();
    void buttonUpdate(int);
    void accept();

signals:
    void loadObject(std::string, double, double, double, double, double, double,double);


protected:
    std::string fileName;
    std::vector< std::string > *list_object;

};

} // namespace qt

} // namespace gui

} // namespace sofa

#endif
