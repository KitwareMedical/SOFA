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
#ifndef SOFA_TUTORIALSELECTOR_H
#define SOFA_TUTORIALSELECTOR_H

#ifdef SOFA_QT4
#include <Q3ListView>
#include <Q3ListViewItem>
typedef Q3ListViewItem QListViewItem;
#include <QKeyEvent>
#else
#include <qlistview.h>
#include <qevent.h>
typedef QListView Q3ListView;
typedef QListViewItem Q3ListViewItem;
#endif

#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/SetDirectory.h>


//Tinyxml library
#include <tinyxml.h>
#include <tinystr.h>

#include <map>

namespace sofa
{

namespace gui
{

namespace qt
{



class TutorialSelector : public Q3ListView
{


    struct Category
    {
        Category() {};
        Category(const std::string &n, const std::string &xml, const std::string &html)
            :name(n), xmlFilename(xml), htmlFilename(html) {};
        std::string name;
        std::string xmlFilename;
        std::string htmlFilename;
        bool operator<(const Category& c) const  { return (name < c.name); }
        bool operator!=(const Category& c) const  { return (name != c.name); }
        bool operator==(const Category& c) const  { return (name == c.name); }
    };

    struct Tutorial
    {
        Tutorial() {};
        Tutorial(const std::string &n, const std::string &scene, const std::string &html)
            :name(n), sceneFilename(scene), htmlFilename(html)
        {
            if ( !sceneFilename.empty() && sofa::helper::system::DataRepository.findFile (sceneFilename) )
                sceneFilename = sofa::helper::system::DataRepository.getFile ( sceneFilename );
            if ( !htmlFilename.empty()  && sofa::helper::system::DataRepository.findFile (htmlFilename) )
                htmlFilename = sofa::helper::system::DataRepository.getFile ( htmlFilename );
        };
        std::string name;
        std::string sceneFilename;
        std::string htmlFilename;
        bool operator<(const Tutorial& t) const  { return (name < t.name); }
        bool operator!=(const Tutorial& t) const  { return (name != t.name); }
        bool operator==(const Tutorial& t) const  { return (name == t.name); }
    };

    Q_OBJECT
public:
    TutorialSelector( QWidget* parent = 0);
    void init();


    void keyPressEvent ( QKeyEvent * e );
    void usingScene(const std::string &filename);
    std::list< std::string > getCategories() const;
public  slots:
    void openCategory(const QString&);

#ifdef SOFA_QT4
    void changeRequested( Q3ListViewItem * );
#else
    void changeRequested( QListViewItem * );
#endif
signals:
    void openCategory(const std::string &name);
    void openTutorial(const std::string &filename);
    void openHTML(const std::string &filename);

protected:

    void openCategory(const Category&);
    void openTutorial(const Tutorial&);

    void loadTutorials(const std::string &fileTutorials);
    void openNode(TiXmlNode* node, Q3ListViewItem *parent=NULL, bool isRoot=false);
    void openAttribute(TiXmlElement* element,  Q3ListViewItem *item);

    std::map< Q3ListViewItem *, Category> itemToCategory;
    std::map< Q3ListViewItem *, Tutorial> itemToTutorial;
    std::multimap<Category, Tutorial> listTutoFromFile;

    Category currentCategory;
};

}
}
}

#endif
