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
#include "SofaMouseManager.h"
#include "ui_MouseManager.h"
#include "QMouseOperations.h"

#include <sofa/gui/MouseOperations.h>
#include <sofa/gui/OperationFactory.h>

#include <iostream>
#ifndef SOFA_QT4
#include <qlineedit.h>
#include <qcombobox.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qlayout.h>
#endif

namespace sofa
{
namespace gui
{
namespace qt
{
SofaMouseManager::SofaMouseManager()
    : gui(new Ui_MouseManager)
{
    gui->setupUi(this);

    connect( gui->LeftOperationCombo,   SIGNAL(activated(int)), this, SLOT( selectOperation(int) ));
    connect( gui->MiddleOperationCombo, SIGNAL(activated(int)), this, SLOT( selectOperation(int) ));
    connect( gui->RightOperationCombo,  SIGNAL(activated(int)), this, SLOT( selectOperation(int) ));

    RegisterOperation("Attach").add< QAttachOperation >();
    RegisterOperation("AddFrame").add< AddFrameOperation >();
    RegisterOperation("SaveCameraViewPoint").add< QAddRecordedCameraOperation >();
    RegisterOperation("StartNavigation").add< QStartNavigationOperation >();
    RegisterOperation("Fix")   .add< QFixOperation  >();
    RegisterOperation("Incise").add< QInciseOperation  >();
    RegisterOperation("Remove").add< QTopologyOperation  >();
    RegisterOperation("Suture").add< QAddSutureOperation >();
    RegisterOperation("ConstraintAttach").add< ConstraintAttachOperation >();
}

SofaMouseManager::~SofaMouseManager()
{
}

void SofaMouseManager::updateContent()
{
    gui->LeftOperationCombo->clear();
    gui->MiddleOperationCombo->clear();
    gui->RightOperationCombo->clear();
    mapIndexOperation.clear();

    if (mapIndexOperation.empty())
    {
        const OperationFactory::RegisterStorage &registry = OperationFactory::getInstance()->registry;

        int idx=0;
        for (OperationFactory::RegisterStorage::const_iterator it=registry.begin(); it!=registry.end(); ++it)
        {
            gui->LeftOperationCombo  ->insertItem(QString(OperationFactory::GetDescription(it->first).c_str()));
            gui->MiddleOperationCombo->insertItem(QString(OperationFactory::GetDescription(it->first).c_str()));
            gui->RightOperationCombo ->insertItem(QString(OperationFactory::GetDescription(it->first).c_str()));

            if (OperationFactory::GetDescription(it->first) == OperationFactory::GetDescription(usedOperations[LEFT]))
                gui->LeftOperationCombo->setCurrentItem(idx);
            if (OperationFactory::GetDescription(it->first) == OperationFactory::GetDescription(usedOperations[MIDDLE]))
                gui->MiddleOperationCombo->setCurrentItem(idx);
            if (OperationFactory::GetDescription(it->first) == OperationFactory::GetDescription(usedOperations[RIGHT]))
                gui->RightOperationCombo->setCurrentItem(idx);

            mapIndexOperation.insert(std::make_pair(idx++, it->first));
        }
    }
}

void SofaMouseManager::setPickHandler(PickHandler *picker)
{
    pickHandler=picker;
    updateContent();
    updateOperation(LEFT,   "Attach");
    updateOperation(MIDDLE, "Incise");
    updateOperation(RIGHT,  "Remove");
}


void SofaMouseManager::selectOperation(int operation)
{
    QComboBox *combo = (QComboBox*)(sender());
    const std::string operationName=mapIndexOperation[operation];

    if      (combo == gui->LeftOperationCombo)   updateOperation(LEFT,   operationName);
    else if (combo == gui->MiddleOperationCombo) updateOperation(MIDDLE, operationName);
    else if (combo == gui->RightOperationCombo)  updateOperation(RIGHT,  operationName);
}

void SofaMouseManager::updateOperation(  sofa::component::configurationsetting::MouseButtonSetting* setting)
{
    //By changing the operation, we delete the previous operation
    Operation *operation=pickHandler->changeOperation( setting);
    updateOperation(operation);
}

void SofaMouseManager::updateOperation( MOUSE_BUTTON button, const std::string &id)
{
    //By changing the operation, we delete the previous operation
    Operation *operation=pickHandler->changeOperation( button, id);
    updateOperation(operation);
}


void SofaMouseManager::updateOperation( Operation* operation)
{
    if (!operation || operation->getMouseButton()==NONE ) return;
    usedOperations[operation->getMouseButton()] = operation->getId();

    QWidget* qoperation=dynamic_cast<QWidget*>(operation);
    if (!qoperation) return;

    switch(operation->getMouseButton())
    {
    case LEFT:
    {
#ifdef SOFA_QT4
        gui->LeftButton->layout()->addWidget(qoperation);
#else
        gui->LeftButton->layout()->add(qoperation);
#endif
        break;
    }
    case MIDDLE:
    {
#ifdef SOFA_QT4
        gui->MiddleButton->layout()->addWidget(qoperation);
#else
        gui->MiddleButton->layout()->add(qoperation);
#endif
        break;
    }
    case RIGHT:
    {
#ifdef SOFA_QT4
        gui->RightButton->layout()->addWidget(qoperation);
#else
        gui->RightButton->layout()->add(qoperation);
#endif
        break;
    }
    default:
    {
    }
    }


}

}
}
}

