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
#ifndef SOFA_GLOBALMODIFICATION_H
#define SOFA_GLOBALMODIFICATION_H

#include <sofa/core/objectmodel/Base.h>

#include "ModifierCondition.h"
#include "GraphHistoryManager.h"

#ifdef SOFA_QT4
#include <Q3Header>
#include <QComboBox>
#include <QLineEdit>
#include <QStringList>
#include <QCheckBox>
#else
#include <qstringlist.h>
#include <qheader.h>
#include <qlineedit.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#endif


namespace sofa
{

namespace gui
{

namespace qt
{




/// Qt Widget applying modifying the value of a given data among a selection of component
class GlobalModification : public QWidget
{
    Q_OBJECT
public:
    typedef helper::vector< sofa::core::objectmodel::Base* > InternalStorage;
    GlobalModification(const InternalStorage &c, GraphHistoryManager* historyManager);
    ~GlobalModification();
public slots:
    void applyGlobalModification();
    void useAliases(bool);
signals:
    void displayMessage(const std::string &message);
protected:
    QCheckBox *aliasEnable;
    QComboBox *dataNameSelector;
    QLineEdit *valueModifier;

    QStringList listDataName;
    QStringList listDataAliases;


    InternalStorage components;
    GraphHistoryManager* historyManager;
    helper::vector< ModifierCondition* > conditions;
};

}
}
}

#endif
