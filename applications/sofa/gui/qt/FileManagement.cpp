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

#include "FileManagement.h"
#include <iostream>

namespace sofa
{
namespace gui
{
namespace qt
{


#ifndef SOFA_QT4
typedef QFileDialog Q3FileDialog;
#include <qdir.h>
#include <qfileinfo.h>
#else
#include <QDir>
#endif

QString getExistingDirectory ( QWidget* parent, const QString & dir, const char * name, const QString & caption)
{
#ifdef SOFA_QT4
    QFileDialog::Options options = QFileDialog::ShowDirsOnly;
    //	options |= QFileDialog::DontUseNativeDialog;
    options |= QFileDialog::DontUseSheet;
    return QFileDialog::getExistingDirectory ( parent, name?QString(name):caption, dir, options );
#else
    return Q3FileDialog::getExistingDirectory( dir, parent, name, caption );
#endif
};

QString getOpenFileName ( QWidget* parent, const QString & startWith, const QString & filter, const char * name, const QString & caption, QString * selectedFilter )
{
#ifdef SOFA_QT4
    QFileDialog::Options options = 0;
    //	options |= QFileDialog::DontUseNativeDialog;
    options |= QFileDialog::DontUseSheet;
    return QFileDialog::getOpenFileName ( parent, name?QString(name):caption, startWith, filter, selectedFilter, options );
#else
    return Q3FileDialog::getOpenFileName ( startWith, filter, parent, name, caption, selectedFilter );
#endif
};

QString getSaveFileName ( QWidget* parent, const QString & startWith, const QString & filter, const char * name, const QString & caption, QString * selectedFilter )
{
#ifdef SOFA_QT4
    QFileDialog::Options options = 0;
    //	options |= QFileDialog::DontUseNativeDialog;
    options |= QFileDialog::DontUseSheet;
    return QFileDialog::getSaveFileName ( parent, name?QString(name):caption, startWith, filter, selectedFilter, options );
#else
    return Q3FileDialog::getSaveFileName ( startWith, filter, parent, name, caption, selectedFilter );
#endif
};

void getFilesInDirectory( const QString &p, std::vector< QString > &files, bool recursive, const std::vector< QString > &filter )
{
    QString path=p;
    if (path.endsWith("/"))
    {
        int slash=path.find('/',-1);
        path.truncate(slash);
    }
    else if (path.endsWith("\\"))
    {
        int slash=path.find('\\',-1);
        path.truncate(slash);
    }

    QDir d(path);

    d.setFilter( QDir::Dirs | QDir::Hidden | QDir::NoSymLinks );

    std::vector< QString > subDir;

    const QFileInfoList &listDirectories =
#ifdef SOFA_QT4
        d.entryInfoList();
    QStringList filters;
    for (unsigned int i=0; i<filter.size(); ++i)
        filters << filter[i];

    d.setNameFilters(filters);
    for (int j = 0; j < listDirectories.size(); ++j)
    {
        QFileInfo fileInfo=listDirectories.at(j);
#else
        *(d.entryInfoList());
    QString filters;
    for (unsigned int i=0; i<filter.size(); ++i)
        filters+=  filter[i] + QString(" ");

    d.setNameFilter(filters);
    QFileInfoListIterator itDir( listDirectories );
    while ( (itDir.current()) != 0 )
    {

        QFileInfo fileInfo=*(itDir.current());
#endif
        subDir.push_back(fileInfo.fileName());
#ifndef SOFA_QT4
        ++itDir;
#endif
    }

    d.setFilter( QDir::Files | QDir::Hidden | QDir::NoSymLinks );

    const QFileInfoList &listFiles =
#ifdef SOFA_QT4
        d.entryInfoList();
    for (int j = 0; j < listFiles.size(); ++j)
    {
        QFileInfo fileInfo=listFiles.at(j);
#else
        *(d.entryInfoList());
    QFileInfoListIterator itFile( listFiles );
    while ( (itFile.current()) != 0 )
    {
        QFileInfo fileInfo=*(itFile.current());
#endif
        files.push_back(path+QString("/")+fileInfo.fileName());
#ifndef SOFA_QT4
        ++itFile;
#endif
    }

    if (recursive)
    {
        for (unsigned int i=0; i<subDir.size(); ++i)
        {
            if (subDir[i].left(1) == QString(".")) continue;
            if (subDir[i] == QString("OBJ"))       continue;

            QString nextDir=path+QString("/")+subDir[i];
            getFilesInDirectory(nextDir, files, recursive, filter);
        }
    }
}

} // namespace qt

} // namespace gui

} // namespace sofa

