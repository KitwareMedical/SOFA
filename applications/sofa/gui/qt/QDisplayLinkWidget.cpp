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
#include "QDisplayLinkWidget.h"
#include "ModifyObject.h"


#ifdef SOFA_QT4
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <Q3GroupBox>
#include <QLabel>
#include <QValidator>
#else
#include <qlayout.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qvalidator.h>
#endif

#define TEXTSIZE_THRESHOLD 45

namespace sofa
{

using namespace core::objectmodel;
using namespace sofa::component::misc;

namespace gui
{

namespace qt
{

QDisplayLinkWidget::QDisplayLinkWidget(QWidget* parent,
        BaseLink* link,
        const ModifyObjectFlags& flags)
    : Q3GroupBox(parent),
      link_(link),
      linkinfowidget_(NULL),
      linkwidget_(NULL),
      numWidgets_(0)
{
    if(link_ == NULL)
    {
        return;
    }

    const std::string label_text = link_->getHelp();

    if (label_text != "TODO")
    {
        linkinfowidget_ = new QDisplayLinkInfoWidget(this,label_text,link_,flags.LINKPATH_MODIFIABLE_FLAG);
		linkinfowidget_->setContentsMargins(0, 0, 0, 0);
		linkinfowidget_->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        
		numWidgets_ += linkinfowidget_->getNumLines()/3;
    }

	setToolTip(link->getHelp());

    LinkWidget::CreatorArgument dwarg;
    dwarg.name =  link_->getName();
    dwarg.link = link_;
    dwarg.parent = this;
    dwarg.readOnly = (!link_->storePath() && flags.READONLY_FLAG);

    linkwidget_= LinkWidget::CreateLinkWidget(dwarg);

    if (linkwidget_ == 0)
    {
        linkwidget_ = new QLinkSimpleEdit(this,dwarg.link->getName().c_str(), dwarg.link);
        linkwidget_->createWidgets();
        linkwidget_->setEnabled( !(dwarg.readOnly) );
        assert(linkwidget_ != NULL);
    }

	if(linkwidget_->layout())
	{
		linkwidget_->layout()->setAlignment(Qt::AlignCenter);
		linkwidget_->layout()->setContentsMargins(2, 2, 2, 2);
	}

	linkwidget_->setContentsMargins(0, 0, 0, 0);
	linkwidget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);

	const std::string valuetype = link_->getValueTypeString();
    if (!valuetype.empty())
        linkwidget_->setToolTip(valuetype.c_str());

	//std::cout << "WIDGET created for link " << dwarg.link << " : " << dwarg.name << " : " << dwarg.link->getValueTypeString() << std::endl;
	numWidgets_ += linkwidget_->sizeWidget();
	connect(linkwidget_,SIGNAL(WidgetDirty(bool)), this, SIGNAL ( WidgetDirty(bool) ) );
	connect(this, SIGNAL( WidgetUpdate() ), linkwidget_, SLOT( updateWidgetValue() ) );
	connect(this, SIGNAL( LinkUpdate() ), linkwidget_, SLOT(updateLinkValue() ) );
	connect(linkwidget_,SIGNAL(LinkOwnerDirty(bool)),this,SIGNAL(LinkOwnerDirty(bool)) );

	if(flags.PROPERTY_WIDGET_FLAG)
    {
		QWidget* refreshWidget = new QWidget(this);
		refreshWidget->setFixedSize(QSize(16, 16));
        QPushButton *refresh = new QPushButton(RefreshIcon(), "", refreshWidget);
        refresh->setHidden(true);
        refresh->setFixedSize(QSize(16, 16));
		refresh->setContentsMargins(0, 0, 0, 0);

        ++numWidgets_;

        connect(linkwidget_,SIGNAL(WidgetDirty(bool)), refresh, SLOT ( setVisible(bool) ) );
        connect(refresh, SIGNAL(clicked()), this, SLOT(UpdateLink()));
        connect(refresh, SIGNAL(clicked(bool)), refresh, SLOT(setVisible(bool)));

		setStyleSheet("QGroupBox{border:0;}");
        setContentsMargins(0, 0, 0, 0);
		setInsideMargin(0);
        setInsideSpacing(0);

        setColumns(numWidgets_);
    }
    else
	{
		setTitle(link_->getName().c_str());
		setInsideMargin(4);
		setInsideSpacing(2);

		setColumns(numWidgets_); //linkwidget_->numColumnWidget()
	}
}

void QDisplayLinkWidget::UpdateLink()
{
    emit LinkUpdate();
}

void QDisplayLinkWidget::UpdateWidgets()
{
    emit WidgetUpdate();
}

QLinkSimpleEdit::QLinkSimpleEdit(QWidget* parent, const char* name, BaseLink* link)
    : LinkWidget(parent,name,link)
{
}

bool QLinkSimpleEdit::createWidgets()
{
    QString str  = QString( getBaseLink()->getValueString().c_str() );
    QLayout* layout = new QHBoxLayout(this);
    if( str.length() > TEXTSIZE_THRESHOLD )
    {
        innerWidget_.type = TEXTEDIT;
        innerWidget_.widget.textEdit = new QTextEdit(this);
		innerWidget_.widget.textEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
		innerWidget_.widget.textEdit->setContentsMargins(0, 0, 0, 0);
		innerWidget_.widget.textEdit->setText(str);
        connect(innerWidget_.widget.textEdit , SIGNAL( textChanged() ), this, SLOT ( update() ) );
        layout->add(innerWidget_.widget.textEdit);
    }
    else
    {
        innerWidget_.type = LINEEDIT;
        innerWidget_.widget.lineEdit  = new QLineEdit(this);
		innerWidget_.widget.lineEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
		innerWidget_.widget.lineEdit->setContentsMargins(0, 0, 0, 0);
        innerWidget_.widget.lineEdit->setText(str);
        connect( innerWidget_.widget.lineEdit, SIGNAL(textChanged(const QString&)), this, SLOT( update() ) );
        layout->add(innerWidget_.widget.lineEdit);
    }
	layout->setContentsMargins(0, 0, 0, 0);
	layout->setSpacing(0);

    return true;
}

void QLinkSimpleEdit::readFromLink()
{
    QString str = QString( getBaseLink()->getValueString().c_str() );
    if(innerWidget_.type == TEXTEDIT)
    {
        innerWidget_.widget.textEdit->setText(str);
    }
    else if(innerWidget_.type == LINEEDIT)
    {
        innerWidget_.widget.lineEdit->setText(str);
    }
}

void QLinkSimpleEdit::writeToLink()
{
    if(getBaseLink())
    {
        std::string value;
        if( innerWidget_.type == TEXTEDIT)
        {
            value = innerWidget_.widget.textEdit->text().ascii();
        }
        else if( innerWidget_.type == LINEEDIT)
        {
            value = innerWidget_.widget.lineEdit->text().ascii();
        }
        getBaseLink()->read(value);
    }
}

} // namespace qt

} // namespace gui

} // namespace sofa
