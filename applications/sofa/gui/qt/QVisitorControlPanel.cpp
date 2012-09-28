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
#include "QVisitorControlPanel.h"
#include <sofa/simulation/common/Visitor.h>

#ifdef SOFA_QT4
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#else
#include <qpushbutton.h>
#include <qcheckbox.h>
#include <qlabel.h>
#include <qlayout.h>
#endif

#if !defined(INFINITY)
#define INFINITY 9.0e10
#endif


namespace sofa
{
namespace gui
{
namespace qt
{

QVisitorControlPanel::QVisitorControlPanel(QWidget* parent): QWidget(parent)
{
    QVBoxLayout *vbox=new QVBoxLayout(this,0,0);


    //Filter the results to quickly find a visitor
    QWidget *filterResult = new QWidget(this);
    QHBoxLayout *hboxFilers=new QHBoxLayout(filterResult);

    textFilter = new QLineEdit(filterResult);
    QPushButton *findFilter = new QPushButton(QString("Find"), filterResult);
    findFilter->setAutoDefault(true);

    hboxFilers->addWidget(new QLabel(QString("Focus on:"), filterResult));
    hboxFilers->addWidget(textFilter);
    hboxFilers->addWidget(findFilter);

    connect(findFilter, SIGNAL(clicked()), this, SLOT(filterResults()));
    connect(textFilter, SIGNAL(returnPressed()), this, SLOT(filterResults()));

    //Clear Button
    QPushButton *clearButton = new QPushButton(QString("Clear"), this);
    clearButton->setAutoDefault(false);
    connect(clearButton, SIGNAL(clicked()), this, SIGNAL(clearGraph()));


    //Parameters to configure the export of the state vectors
    QWidget *exportStateParameters = new QWidget(this);

    QHBoxLayout *hboxParameters=new QHBoxLayout(exportStateParameters);

    QCheckBox *activation=new QCheckBox(QString("Trace State Vector"), exportStateParameters);

    spinIndex = new WDoubleLineEdit(exportStateParameters, "index");
    spinIndex->setMinValue( (double)-INFINITY );
    spinIndex->setMaxValue( (double)INFINITY );
    spinIndex->setIntValue(sofa::simulation::Visitor::GetFirstIndexStateVector());
    spinIndex->setMaximumWidth(50);
    spinRange = new WDoubleLineEdit(exportStateParameters, "range");
    spinRange->setMinValue( (double)-INFINITY );
    spinRange->setMaxValue( (double)INFINITY );
    spinRange->setIntValue(sofa::simulation::Visitor::GetRangeStateVector());
    spinRange->setMaximumWidth(50);

    connect(activation, SIGNAL(toggled(bool)), this, SLOT(activateTraceStateVectors(bool)));
    connect(spinIndex, SIGNAL(lostFocus()), this, SLOT(changeFirstIndex()));
    connect(spinRange, SIGNAL(lostFocus()), this, SLOT(changeRange()));


    hboxParameters->addWidget(activation);
    hboxParameters->addItem(new QSpacerItem(100,10));
    hboxParameters->addWidget(new QLabel(QString("First Index"), exportStateParameters));
    hboxParameters->addWidget(spinIndex);

    hboxParameters->addWidget(new QLabel(QString("Range"), exportStateParameters));
    hboxParameters->addWidget(spinRange);

    hboxParameters->addStretch();

    //Configure the main window
    vbox->addWidget(filterResult);
    vbox->addWidget(exportStateParameters);
    vbox->addWidget(clearButton);

    activateTraceStateVectors(sofa::simulation::Visitor::IsExportStateVectorEnabled());

}

void QVisitorControlPanel::activateTraceStateVectors(bool active)
{
    sofa::simulation::Visitor::EnableExportStateVector(active);
    spinIndex->setEnabled(active);
    spinRange->setEnabled(active);
}
void QVisitorControlPanel::changeFirstIndex()
{
    WDoubleLineEdit *w=(WDoubleLineEdit *) sender();
    int value=w->getIntValue();
    changeFirstIndex(value);
}
void QVisitorControlPanel::changeRange()
{
    WDoubleLineEdit *w=(WDoubleLineEdit *) sender();
    int value=w->getIntValue();
    changeRange(value);
}

void QVisitorControlPanel::changeFirstIndex(int idx)
{
    sofa::simulation::Visitor::SetFirstIndexStateVector(idx);
}
void QVisitorControlPanel::changeRange(int range)
{
    sofa::simulation::Visitor::SetRangeStateVector(range);
}
void QVisitorControlPanel::filterResults()
{
    emit focusOn(textFilter->text());
}



} // qt
} // gui
} //sofa


