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
#ifndef SOFA_GUI_VIEWER_QT_QMENUFILESRECENTLYOPENED_H
#define SOFA_GUI_VIEWER_QT_QMENUFILESRECENTLYOPENED_H

#include <sofa/gui/FilesRecentlyOpenedManager.h>
#include "SofaGUIQt.h"
#ifdef SOFA_QT4
#include <QMenu>
#else
#include <qpopupmenu.h>
typedef QPopupMenu QMenu;
#endif
namespace sofa
{
namespace gui
{
namespace qt
{

class SOFA_SOFAGUIQT_API QMenuFilesRecentlyOpened: public FilesRecentlyOpenedManager
{
public:
    QMenuFilesRecentlyOpened(const std::string &configFile):FilesRecentlyOpenedManager(configFile),menuRecentlyOpenedFiles(0) {};
    virtual ~QMenuFilesRecentlyOpened() {if (menuRecentlyOpenedFiles) delete menuRecentlyOpenedFiles;};
    void openFile(const std::string &file);

    QMenu *createWidget(QWidget *parent, const std::string& =std::string("Recently Opened Files ..."));
    QMenu *getMenu() {return menuRecentlyOpenedFiles;};

protected:
    void updateWidget();
    QMenu *menuRecentlyOpenedFiles;


};


}
}
}

#endif
