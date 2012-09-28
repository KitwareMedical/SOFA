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
#include "PieWidget.h"
#include <iostream>

#ifdef SOFA_QT4
#include <QGridLayout>
#include <QStringList>
#include <QHeaderView>
#include <QSplitter>
#include <Qt>
#else
#include <qheader.h>
#include <qlayout.h>
#include <qsplitter.h>
#endif



namespace sofa
{

namespace gui
{

namespace qt
{

std::vector< defaulttype::Vec<3,int> > PieWidget::colorArray;

defaulttype::Vec<3,int> PieWidget::getColor(int i)
{
    defaulttype::Vec<3,int> res=PieWidget::colorArray[i%PieWidget::colorArray.size()];
    float factor=1.0/(1.0+(0.3*(i/PieWidget::colorArray.size())));
    res[0] = (int)(res[0]*factor);
    res[1] = (int)(res[1]*factor);
    res[2] = (int)(res[2]*factor);
    return res;
}

PieWidget::PieWidget(QWidget *parent): QWidget(parent)
{
    if (PieWidget::colorArray.empty())
    {
        colorArray.push_back(  defaulttype::Vec<3,int>(250,125,70) );
        colorArray.push_back(  defaulttype::Vec<3,int>(120,220,110) );
        colorArray.push_back(  defaulttype::Vec<3,int>(215,90,215) );
        colorArray.push_back(  defaulttype::Vec<3,int>(255,210,40) );
        colorArray.push_back(  defaulttype::Vec<3,int>(75,210,210) );
    }
}
void PieWidget::paintEvent( QPaintEvent* )
{
    sizePie = (int)(std::min(this->width(),this->height())*0.95);
    if (data.empty()) return;

    QPainter p( this );

    int initDraw=0;
    int startAngle=0;

    p.setBrush(Qt::SolidPattern);

    for (unsigned int i=0; i<data.size() && i<selection; ++i)
    {
        defaulttype::Vec<3,int> c=PieWidget::getColor(i);
        QColor color(c[0],c[1],c[2]);
        p.setBrush(color);
        int spanAngle=(int)(16*360*data[i].time/totalTime);
        p.drawPie(initDraw,initDraw,sizePie,sizePie,startAngle,spanAngle);
        startAngle+= spanAngle;
    }
}

void PieWidget::setChart( std::vector< dataTime >& value, unsigned int s)
{
    data=value;
    selection=s;
    totalTime=0;
    for (unsigned int i=0; i<value.size() && i<selection; ++i)
    {
        totalTime +=data[i].time;
    }
}

void PieWidget::clear()
{
    data.clear();
    repaint();
}

ChartsWidget::ChartsWidget(const std::string &name, QWidget *parent): QWidget(parent)
{
    QSplitter *splitter=new QSplitter(this);
    splitter->setOrientation(Qt::Horizontal);
    QGridLayout *grid = new QGridLayout(this);
    pie = new PieWidget(splitter);
#ifdef SOFA_QT4
    table = new QTableWidget(0,3,splitter);
    table->horizontalHeader()->setResizeMode(0,QHeaderView::Fixed);
    table->horizontalHeader()->resizeSection(0,30);
    table->horizontalHeader()->setResizeMode(1,QHeaderView::ResizeToContents);
    table->horizontalHeader()->setResizeMode(2,QHeaderView::ResizeToContents);
    QStringList list; list<<"Id" << name.c_str() << "Time";
    table->setHorizontalHeaderLabels(list);
#else
    table = new QTableWidget(0,2,splitter);
    table->horizontalHeader()->setLabel(0,QString(name.c_str()));
    table->horizontalHeader()->setStretchEnabled(true,0);
    table->horizontalHeader()->setLabel(1,QString("Time"));
    table->horizontalHeader()->setStretchEnabled(true,1);
#endif

    grid->addWidget(splitter,0,0);

}


void ChartsWidget::clear()
{
#ifdef SOFA_QT4
    int rows=table->rowCount();
#else
    int rows=table->numRows();
#endif
    for (int i=0; i<rows; ++i) table->removeRow(0);
    pie->clear();
}

void ChartsWidget::setChart( std::vector< dataTime >& value, unsigned int s)
{
    clear();
    pie->setChart(value,s);
    selection=s;
    for (unsigned int i=0; i<value.size() && i<selection; ++i)
    {
#ifdef SOFA_QT4
        table->insertRow(i);
#else
        table->insertRows(i);
#endif
        defaulttype::Vec<3,int> c=PieWidget::getColor(i);
        QColor color(c[0],c[1],c[2]);

        QString text(value[i].name.c_str());
        QString time= QString::number(value[i].time);
        time += QString(" ms");
        if (!value[i].type.empty())
        {
            text+="(";
            text+= QString(value[i].type.c_str());
            text+=")";
        }

#ifdef SOFA_QT4
        QTableWidgetItem *itemColor = new QTableWidgetItem();
        itemColor->setBackgroundColor(color);
        QTableWidgetItem *item = new QTableWidgetItem();
        QTableWidgetItem *itemTime = new QTableWidgetItem();
        table->setItem(i,0, itemColor);
        item->setText(text);
        table->setItem(i,1, item);
        itemTime->setText(time);
        table->setItem(i,2, itemTime);
        table->resizeColumnToContents(1);
        itemColor->setFlags(0);
        item->setFlags(0);
        itemTime->setFlags(0);
#else
        QTableWidgetItem *item = new QTableWidgetItem(table, QTableItem::Never);
        QPixmap p(20,20); p.fill(color);
        item->setPixmap(p);
        QTableWidgetItem *itemTime = new QTableWidgetItem(table, QTableItem::Never);
        item->setText(text);
        table->setItem(i,0, item);
        itemTime->setText(time);
        table->setItem(i,1, itemTime);
        table->adjustColumn(0);
#endif
    }
    pie->repaint();

}

}
}
}
