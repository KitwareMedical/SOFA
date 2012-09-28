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
#ifndef SOFA_VIEWER_H
#define SOFA_VIEWER_H

#include "../../BaseViewer.h"
#include "VisualModelPolicy.h"
#include "../PickHandlerCallBacks.h"
#include "../SofaGUIQt.h"
#include "../SofaVideoRecorderManager.h"

#include <qstring.h>
#include <qwidget.h>

#ifdef SOFA_QT4
#include <QEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QTabWidget>
#include <QTimer>
#else
#include <qevent.h>
#include <qtabwidget.h>
#include <qtimer.h>
#endif

#include <qcursor.h>

namespace sofa
{

namespace gui
{

namespace qt
{

namespace viewer
{

using namespace sofa::defaulttype;
enum
{
    BTLEFT_MODE = 101, BTRIGHT_MODE = 102, BTMIDDLE_MODE = 103,
};




class SOFA_SOFAGUIQT_API SofaViewer : public sofa::gui::BaseViewer
{

public:
    SofaViewer();
    virtual ~SofaViewer();

    /// Optional QTabWidget GUI for a concreate viewer.
    virtual void removeViewerTab(QTabWidget *) {}
    /// Optional QTabWidget GUI for a concreate viewer.
    virtual void configureViewerTab(QTabWidget *) {}

    virtual QWidget* getQWidget()=0;
    virtual QString helpString()=0;

    //*************************************************************
    // QT
    //*************************************************************
    //SLOTS
    virtual void captureEvent();

    // ---------------------- Here are the Keyboard controls   ----------------------
    virtual void keyPressEvent(QKeyEvent * e);
    virtual void keyReleaseEvent(QKeyEvent * e);
    bool isControlPressed() const;
    // ---------------------- Here are the Mouse controls   ----------------------
    virtual void wheelEvent(QWheelEvent *e);
    virtual void mouseMoveEvent ( QMouseEvent *e );
    virtual void mousePressEvent ( QMouseEvent * e);
    virtual void mouseReleaseEvent ( QMouseEvent * e);
    virtual bool mouseEvent(QMouseEvent *e);

protected:
    virtual void redraw();

    QTimer captureTimer;

    bool m_isControlPressed;

    ColourPickingRenderCallBack colourPickingRenderCallBack;

signals:
    virtual void redrawn() = 0;
    virtual void resizeW(int) = 0;
    virtual void resizeH(int) = 0;
};

template < typename VisualModelPolicyType >
class CustomPolicySofaViewer : public VisualModelPolicyType, public sofa::gui::qt::viewer::SofaViewer
{
public:
    using VisualModelPolicyType::load;
    using VisualModelPolicyType::unload;
    CustomPolicySofaViewer() { load(); }
    virtual ~CustomPolicySofaViewer() { unload(); }
protected:
};

typedef CustomPolicySofaViewer< OglModelPolicy > OglModelSofaViewer;


}
}
}
}

#endif
