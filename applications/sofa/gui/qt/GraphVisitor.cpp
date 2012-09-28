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
#include "GraphVisitor.h"

#include <sstream>

#include <tinyxml.cpp>
#include <tinyxmlerror.cpp>
#include <tinystr.cpp>
#include <tinyxmlparser.cpp>

#include <algorithm>

#include <sofa/helper/system/thread/CTime.h>

#ifndef SOFA_QT4
#include <qsplitter.h>
#include <qvaluelist.h>
typedef QListViewItem Q3ListViewItem;
#else
#include <QSplitter>
#include <QList>
#endif
namespace sofa
{

namespace gui
{

namespace qt
{

typedef sofa::helper::system::thread::CTime CTime;

bool cmpTime(const dataTime &a, const dataTime &b) { return a.time > b.time;};
bool GraphVisitor::load(std::string &file)
{
    //Open it using TinyXML
    TiXmlDocument doc;
    doc.Parse(file.c_str());

    //std::cerr << "GRAPH:"<< std::endl << file << std::endl;

    TiXmlHandle hDoc(&doc);
    TiXmlNode* pElem;
    //Getting the root of the file
    pElem=hDoc.FirstChildElement().Element();

    // should always have a valid root but handle gracefully if it does
    if (!pElem) return false;

    totalTime = getTotalTime(pElem);

    componentsTime.clear();
    visitorsTime.clear();


    openNode( pElem, NULL, NULL);


    std::sort(componentsTime.begin(),componentsTime.end(),cmpTime);
    std::sort(visitorsTime.begin(),visitorsTime.end(),cmpTime);

    std::sort(componentsTimeTotal.begin(),componentsTimeTotal.end(),cmpTime);
    std::sort(visitorsTimeTotal.begin(),visitorsTimeTotal.end(),cmpTime);


    if (totalTimeMax<totalTime)
    {
        totalTimeMax=totalTime;
        componentsTimeMax=componentsTime;
        visitorsTimeMax=visitorsTime;
    }

    //        window->pieChart->setChart(visitorsTime, visitorsTime.size());
    window->setCharts(componentsTime,componentsTimeMax,componentsTimeTotal,
            visitorsTime,visitorsTimeMax,visitorsTimeTotal);

    if (!initSize)
    {
        const int sizeLeft = window->graphView->columnWidth(0)+window->graphView->columnWidth(1)+7;
#ifdef SOFA_QT4
        QList< int > listSize;
#else
        QValueList< int > listSize;
#endif
        listSize << sizeLeft
                << window->statsWidget->width()-(sizeLeft-window->graphView->width());
        window->splitterStats->setSizes(listSize);
        initSize=true;
    }
    return true;
}


void GraphVisitor::openAttribute      ( TiXmlElement* element, Q3ListViewItem* item)
{
    if (!element) return;
    TiXmlAttribute* attribute=element->FirstAttribute();
    while (attribute)
    {
        std::string nameOfAttribute(attribute->Name());
        std::string valueOfAttribute(attribute->Value());
        addInformation(item, nameOfAttribute, valueOfAttribute);
        attribute=attribute->Next();
    }

}



void GraphVisitor::openTime      ( TiXmlNode* node, Q3ListViewItem* item)
{
    TiXmlElement* element=node->ToElement();
    if (!element) return;
    TiXmlAttribute* attribute=element->FirstAttribute();
    double timeSec= getTime(attribute);
    double time = 100.0*timeSec/totalTime;
    std::ostringstream s;
    s.setf(std::ios::fixed, std::ios::floatfield);
    s.precision(3);

    s << time << "%";

    TiXmlNode* parent = node->Parent();
    if (parent)
    {
        std::string nodeType = parent->Value();
        if (nodeType == "Component")
        {
            TiXmlAttribute* attribute=parent->ToElement()->FirstAttribute();
            std::string componentName, componentType, componentPtr;
            while (attribute)
            {
                std::string nameOfAttribute(attribute->Name());
                std::string valueOfAttribute(attribute->Value());
                if (nameOfAttribute=="name")
                    componentName=valueOfAttribute;
                else if (nameOfAttribute=="type")
                    componentType=valueOfAttribute;
                else if (nameOfAttribute=="ptr")
                    componentPtr=valueOfAttribute;
                attribute=attribute->Next();
            }
            if (std::find(visitedNode.begin(), visitedNode.end(), componentName) == visitedNode.end())
            {
                dataTime t(timeSec-timeComponentsBelow.back()
                        , componentType, componentName, componentPtr);
                std::vector< dataTime >::iterator it=std::find(componentsTime.begin(),componentsTime.end(),t);
                if (it != componentsTime.end()) it->time += t.time;
                else componentsTime.push_back(t);


                it=std::find(componentsTimeTotal.begin(),componentsTimeTotal.end(),t);
                if (it != componentsTimeTotal.end()) it->time += t.time;
                else componentsTimeTotal.push_back(t);


                visitedNode.push_back(componentName);
            }
        }
        else
        {
            if (nodeType != "Node" && nodeType != "Input" && nodeType != "Output" && nodeType != "Vector" &&
                std::find(visitedNode.begin(), visitedNode.end(),nodeType) == visitedNode.end())
            {
                dataTime t(timeSec, nodeType);
                std::vector< dataTime >::iterator it=std::find(visitorsTime.begin(),visitorsTime.end(),t);
                if (it != visitorsTime.end()) it->time += timeSec;
                else visitorsTime.push_back(t);

                it=std::find(visitorsTimeTotal.begin(),visitorsTimeTotal.end(),t);
                if (it != visitorsTimeTotal.end()) it->time += t.time;
                else visitorsTimeTotal.push_back(t);


                visitedNode.push_back(nodeType);
            }
        }

        if (nodeType == "Node" || nodeType == "Component" || nodeType.rfind("Visitor") == nodeType.size()-7)
            timeComponentsBelow.back() = timeSec;
    }

    addTime(item,  s.str());
}

double GraphVisitor::getTime(TiXmlAttribute* attribute) const
{
    static double conversion=1.0/(double)CTime::getTicksPerSec();
    std::string valueOfAttribute(attribute->Value());
    double result=1000.0*atof(valueOfAttribute.c_str())*conversion;
    return result;
}

double GraphVisitor::getTotalTime(TiXmlNode* node) const
{

    for ( TiXmlNode* child = node->FirstChild(); child != 0; child = child->NextSibling())
    {
        std::string nameOfNode=child->Value();
        if (nameOfNode == "TotalTime")
        {
            TiXmlAttribute* attribute=child->ToElement()->FirstAttribute();
            double total=getTime(attribute);
            std::ostringstream out; out << total;
            attribute->SetValue(out.str().c_str());
            return total;
        }
    }
    return 1;
}

Q3ListViewItem* GraphVisitor::openNode( TiXmlNode* node, Q3ListViewItem* parent, Q3ListViewItem* elementAbove)
{
    if (!node) return NULL;

    unsigned int sizeVisitedNode=visitedNode.size();
    std::string nameOfNode=node->Value();
    // TinyXml API changed in 2.6.0, ELEMENT was replaced with TINYXML_ELEMENT
    // As the version number is not available as a macro, the most portable was is to
    // replace these constants with checks of the return value of ToElement(), ...
    // -- Jeremie A. 02/07/2011
    //int typeOfNode=node->Type();
    Q3ListViewItem *graphNode=NULL;
    if (node->ToDocument())   // case TiXmlNode::DOCUMENT:
    {
    }
    else if (node->ToElement())     // case TiXmlNode::ELEMENT:
    {
        if (nameOfNode == "Time")
        {
            openTime( node, parent);
        }
        else
        {
            graphNode = addNode(parent, elementAbove, nameOfNode);
            openAttribute( node->ToElement(), graphNode);
        }
    }
    else if (node->ToComment())     // case TiXmlNode::COMMENT:
    {
        graphNode = addComment(parent, elementAbove, nameOfNode);
    }
    else if (node->ToText())     // case TiXmlNode::TEXT:
    {
    }
    else if (node->ToDeclaration())     // case TiXmlNode::DECLARATION:
    {
    }
    else     // default:
    {
    }

    Q3ListViewItem *element=NULL;
    timeComponentsBelow.push_back(0);

    for ( TiXmlNode* child = node->FirstChild(); child != 0; child = child->NextSibling())
    {
        element = openNode( child, graphNode, element);
    }
    double t=timeComponentsBelow.back();

    timeComponentsBelow.pop_back();

    if (!timeComponentsBelow.empty()) timeComponentsBelow.back() += t;

    if (sizeVisitedNode != visitedNode.size()) visitedNode.resize(sizeVisitedNode);
    return graphNode;
}


Q3ListViewItem *GraphVisitor::addNode(Q3ListViewItem *parent, Q3ListViewItem *elementAbove, std::string name)
{
    Q3ListViewItem *item=NULL;
    if (!parent)
    {
        //Add a Root
        item=new Q3ListViewItem(graph, QString(name.c_str()));
        item->setOpen(true);
    }
    else
    {
        //Add a child to a node
        item=new Q3ListViewItem(parent,elementAbove, QString(name.c_str()));
    }
    QPixmap*  icon=WindowVisitor::getPixmap(WindowVisitor::getComponentType(name));
    if (icon) item->setPixmap(0,*icon);
    item->setMultiLinesEnabled(true);
    return item;
}

void GraphVisitor::addTime(Q3ListViewItem *element, std::string info)
{
    if (!element) return;
    element->setText(1, QString( info.c_str()));
}

void GraphVisitor::addInformation(Q3ListViewItem *element, std::string name, std::string info)
{
    if (!element) return;
    if (element->text(0) == QString("Node"))
        element->setText(0, QString(info.c_str()));
    else  if (element->text(0) == QString("Component"))
        element->setText(0, QString(info.c_str()));
    else
    {

        QString nameQt = element->text(2);
        QString infoQt = element->text(3);
        if (!nameQt.isEmpty())
        {
            nameQt += QString("\n");
            infoQt += QString("\n");
        }

        if (element->text(0) == QString("Vector")  && name=="value" && !info.empty())
        {
            std::istringstream ss(info);
            std::ostringstream result;
            unsigned int size; ss >> size;
            while (!ss.eof())
            {
                result << "[";
                for (unsigned int i=0; i<size; ++i)
                {
                    float v; ss >> v;
                    result << v;
                    if (i!=size-1) result << " ";
                }
                result << "]\n";
            }
            info = result.str();
        }

        nameQt += QString(name.c_str());
        infoQt += QString(info.c_str());

        if (name != "ptr")
        {
            element->setText(2, nameQt);
            element->setText(3, infoQt);
        }
    }
}

Q3ListViewItem *GraphVisitor::addComment(Q3ListViewItem *element,Q3ListViewItem *elementAbove,  std::string comment)
{
    if (!element) return NULL;
    Q3ListViewItem *result = new Q3ListViewItem(element, elementAbove,QString(comment.c_str()));
    result->setPixmap(0,*WindowVisitor::getPixmap(WindowVisitor::COMMENT));
    result->setSelectable(false);
    result->setMultiLinesEnabled(true);
    return result;
}


}
}
}
