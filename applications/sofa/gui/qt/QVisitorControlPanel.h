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
#ifndef SOFA_GUI_QT_SOFA_GUI_QT_QVISITORCONTROLPANEL_H
#define SOFA_GUI_QT_SOFA_GUI_QT_QVISITORCONTROLPANEL_H

#include <sofa/simulation/common/Node.h>

#ifdef SOFA_QT4
#include <QWidget>
#else
#include <qwidget.h>
#endif

#include "WDoubleLineEdit.h"


namespace sofa
{
namespace gui
{
namespace qt
{

class QVisitorControlPanel : public QWidget
{
    Q_OBJECT
public:
    QVisitorControlPanel(QWidget* parent);

    void changeFirstIndex(int);
    void changeRange(int);
public slots:
    void activateTraceStateVectors(bool);
    void changeFirstIndex();
    void changeRange();
    void filterResults();
signals:
    void focusOn(QString);
    void clearGraph();
protected:
    QLineEdit *textFilter;
    WDoubleLineEdit *spinIndex;
    WDoubleLineEdit *spinRange;
};


} // qt
} // gui
} //sofa

#endif // SOFA_GUI_QT_QVISITORCONTROLPANEL_H

