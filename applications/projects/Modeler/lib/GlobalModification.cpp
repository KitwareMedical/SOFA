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

#include "GlobalModification.h"

#include <sofa/helper/set.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/objectmodel/BaseObject.h>

#ifdef SOFA_QT4
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QCompleter>
#include <QPushButton>
#include <QSpacerItem>
#include <QGroupBox>
#else
#include <qlayout.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <qsizepolicy.h>
#include <qgroupbox.h>
#endif

namespace sofa
{

namespace gui
{

namespace qt
{


GlobalModification::GlobalModification(const InternalStorage &c, GraphHistoryManager* h): components(c), historyManager(h)
{
    setCaption(QString("Global Modifications"));


    helper::set< std::string > allNames;
    helper::set< std::string > allAliases;

    for (InternalStorage::const_iterator it=components.begin(); it!=components.end(); ++it)
    {
        const core::objectmodel::Base *c=(*it);
        const core::objectmodel::Base::VecData& datas=c->getDataFields();
        for (unsigned int d=0; d<datas.size(); ++d)
        {
            allNames.insert(datas[d]->getName());
        }

        const core::objectmodel::Base::MapData& aliases=c->getDataAliases();
        for (core::objectmodel::Base::MapData::const_iterator it=aliases.begin(); it!=aliases.end(); ++it)
        {
            allAliases.insert(it->first);
        }

    }

    for (helper::set< std::string >::const_iterator it=allNames.begin(); it!=allNames.end(); ++it)
        listDataName << it->c_str();

    for (helper::set< std::string >::const_iterator it=allAliases.begin(); it!=allAliases.end(); ++it)
        listDataAliases << it->c_str();

    //Creation of the GUI
    QVBoxLayout *globalLayout = new QVBoxLayout(this);




    //***********************************************************************************
    //Selection of the Data Name
    QWidget *dataNameWidget = new QWidget(this);
    QHBoxLayout *nameLayout = new QHBoxLayout(dataNameWidget);

    dataNameSelector = new QComboBox();

    dataNameSelector->setEditable(true);
    dataNameSelector->setAutoCompletion(true);

    dataNameSelector->insertStringList(listDataName);
#ifdef SOFA_QT4
    dataNameSelector->setCompleter(new QCompleter(listDataName,dataNameSelector));
#endif

    nameLayout->addWidget(new QLabel(QString("Select the Data to modify"),dataNameWidget));
    nameLayout->addWidget(dataNameSelector);

    //***********************************************************************************
    //New value
    QWidget *dataValueWidget = new QWidget(this);
    QHBoxLayout *valueLayout = new QHBoxLayout(dataValueWidget);

    valueModifier = new QLineEdit(dataValueWidget);


    valueLayout->addWidget(new QLabel(QString("New value: "),dataValueWidget));
    valueLayout->addWidget(valueModifier);



    //***********************************************************************************
    //Button Panel
    QWidget *buttonWidget = new QWidget(this);
    QHBoxLayout *buttonLayout = new QHBoxLayout(buttonWidget);

    QPushButton *modifyButton = new QPushButton(QString("Modify"), buttonWidget);
    QPushButton *cancelButton = new QPushButton(QString("Cancel"), buttonWidget);

    modifyButton->setAutoDefault(true);

    buttonLayout->addItem(new QSpacerItem(1,1,QSizePolicy::Expanding,QSizePolicy::Minimum));
    buttonLayout->addWidget(modifyButton);
    buttonLayout->addWidget(cancelButton);

    //***********************************************************************************
    //Control Panel

    QGroupBox *controlPanel = new QGroupBox(this);
    QVBoxLayout *controlLayout = new QVBoxLayout(controlPanel);

    //--------------------
    //Alias Enable
    QWidget *dataAliasEnableWidget = new QWidget(this);
    QHBoxLayout *aliasLayout = new QHBoxLayout(dataAliasEnableWidget);

    aliasEnable= new QCheckBox(QString("Using also Aliases for Data name"),dataAliasEnableWidget);

    aliasLayout->addWidget(aliasEnable);



    //--------------------
    //Name Condition
    QNamingModifierCondition *nameCondition = new QNamingModifierCondition();
    conditions.push_back(nameCondition);

    //--------------------
    //Value Condition
    QValueModifierCondition *valueCondition = new QValueModifierCondition();
    conditions.push_back(valueCondition);

    //***********************************************************************************
    // Setup the layout
    controlLayout->addWidget(dataAliasEnableWidget);
    controlLayout->addWidget(nameCondition);
    controlLayout->addWidget(valueCondition);


    globalLayout->addWidget(dataNameWidget);
    globalLayout->addWidget(controlPanel);
    globalLayout->addWidget(dataValueWidget);
    globalLayout->addWidget(buttonWidget);


    //***********************************************************************************
    //Do the connections
    connect(aliasEnable, SIGNAL(toggled(bool)), this, SLOT(useAliases(bool)));
    connect(modifyButton, SIGNAL(clicked()), this, SLOT(applyGlobalModification()));
    connect(cancelButton, SIGNAL(clicked()), this, SLOT(close()));
    connect(valueModifier, SIGNAL(returnPressed()), this, SLOT(applyGlobalModification()));
}

GlobalModification::~GlobalModification()
{
    delete dataNameSelector;
    delete valueModifier;
}


void GlobalModification::applyGlobalModification()
{
    QStringList *list;

    if (aliasEnable->isChecked ()) list = &listDataAliases;
    else   list = &listDataName;


    QString name=dataNameSelector->currentText();
    if (list->find(name) == list->end())
    {
        const std::string message="Data " + std::string(name.ascii()) + " does not exist!";
        emit displayMessage(message);
        close();
        return;
    }

    std::string n=name.ascii();
    std::string v=valueModifier->text().ascii();

    for (InternalStorage::const_iterator it=components.begin(); it!=components.end(); ++it)
    {
        core::objectmodel::Base *c=(*it);

        const std::vector< core::objectmodel::BaseData* > &data=c->findGlobalField(n);
        if (!data.empty())
        {

            if (historyManager) historyManager->beginModification(c);


            for (unsigned int i=0; i<data.size(); ++i)
            {
                bool conditionsAccepted=true;
                for (unsigned int cond=0; cond<conditions.size() && conditionsAccepted; ++cond) conditionsAccepted &= conditions[cond]->verify(c,data[i]);
                if (conditionsAccepted) data[i]->read(v);
            }
            if (historyManager) historyManager->endModification(c);
        }
    }

    close();
}


void GlobalModification::useAliases(bool b)
{
    const QString currentText=dataNameSelector->currentText();
    dataNameSelector->clear();

    QStringList *list;

    if (b) list = &listDataAliases;
    else   list = &listDataName;


#ifdef SOFA_QT4
    if (dataNameSelector->completer()) delete dataNameSelector->completer();
    dataNameSelector->setCompleter(new QCompleter(*list,dataNameSelector));
#endif
    dataNameSelector->insertStringList(*list);

    dataNameSelector->setCurrentText(currentText);
}

}
}
}

