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
#ifndef SOFA_GUI_FILEMANAGEMENT_H
#define SOFA_GUI_FILEMANAGEMENT_H

#include "SofaGUIQt.h"
#include <vector>

#ifdef SOFA_QT4
#include <Q3FileDialog>
#include <QFileDialog>
#else
#include <qfiledialog.h>
#endif

namespace sofa
{

namespace gui
{

namespace qt
{


#ifndef SOFA_QT4
typedef QFileDialog Q3FileDialog;
#endif

QString SOFA_SOFAGUIQT_API getExistingDirectory ( QWidget* parent, const QString & dir = QString(), const char * name = 0, const QString & caption = QString() );

QString SOFA_SOFAGUIQT_API getOpenFileName ( QWidget* parent, const QString & startWith = QString(), const QString & filter = QString(), const char * name = 0, const QString & caption = QString(), QString * selectedFilter = 0 );

QString SOFA_SOFAGUIQT_API getSaveFileName ( QWidget* parent, const QString & startWith = QString(), const QString & filter = QString(), const char * name = 0, const QString & caption = QString(), QString * selectedFilter = 0 );


void SOFA_SOFAGUIQT_API getFilesInDirectory( const QString &path, std::vector< QString > &files, bool recursive=true, const std::vector< QString > &filter=std::vector< QString >() );

} // namespace qt

} // namespace gui

} // namespace sofa

#endif // SOFA_GUI_VIEWER_REALGUI_H
