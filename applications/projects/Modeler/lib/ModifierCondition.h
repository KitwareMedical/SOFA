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
#ifndef SOFA_MODIFIERCONDITION_H
#define SOFA_MODIFIERCONDITION_H

#include <sofa/core/objectmodel/Base.h>

#ifdef SOFA_QT4
#include <QCheckBox>
#include <QLineEdit>
#include <QComboBox>
#include <QStringList>
#else
#include <qcheckbox.h>
#include <qlineedit.h>
#include <qcombobox.h>
#include <qstringlist.h>
#endif

namespace sofa
{

namespace gui
{

namespace qt
{


struct ModifierCondition
{
    virtual ~ModifierCondition() {};
    virtual bool verify(core::objectmodel::Base* c, core::objectmodel::BaseData* d) const =0;
    virtual bool isActive() const=0;
};




class QNamingModifierCondition: public QWidget, public ModifierCondition
{
public:
    QNamingModifierCondition(QWidget *parent=0);

    bool verify(core::objectmodel::Base* c, core::objectmodel::BaseData* d) const;
    bool isActive() const {return activated->isChecked();}
    std::string getValue() const {return entryName->text().ascii();}
protected:
    QCheckBox *activated;
    QLineEdit *entryName;
    QComboBox *criteriaSelector;
    QStringList criteriaList;
};


class QValueModifierCondition: public QWidget, public ModifierCondition
{
public:
    QValueModifierCondition(QWidget *parent=0);

    bool verify(core::objectmodel::Base* c, core::objectmodel::BaseData* d) const;

    bool isActive() const {return activated->isChecked();}
    std::string getValue() const {return value->text().ascii();}
protected:
    QCheckBox *activated;
    QLineEdit *value;
    QComboBox *criteriaSelector;
    QStringList criteriaList;
};



}
}
}

#endif
