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
#ifndef SOFA_GUI_QT_MODELVIEWTABLEDATACONTAINER_H
#define SOFA_GUI_QT_MODELVIEWTABLEDATACONTAINER_H

#include "SimpleDataWidget.h"
#include "StructDataWidget.h"
#ifdef TODOTOPO
#include <sofa/component/topology/PointSubsetData.h>
#endif
#include <sofa/component/topology/TopologyData.h>

#include "QModelViewTableUpdater.h"
#include "QDisplayDataWidget.h"

#include <QVBoxLayout>
#include <QTableView>
#include <QStandardItemModel>

#include <sofa/helper/deque.h>

namespace sofa
{

namespace gui
{

namespace qt
{

enum
{
    TYPE_SINGLE = 0,
    TYPE_VECTOR = 1,
    TYPE_STRUCT = 2,
};

template<class T, int TYPE>
class flat_data_trait;

template<class T>
class default_flat_data_trait : public flat_data_trait< T,  (  (struct_data_trait<T>::NVAR >1 ) ?      TYPE_STRUCT :     (       (vector_data_trait<T>::NDIM > 0) ?       TYPE_VECTOR :       TYPE_SINGLE ) ) > {};



template<class T> inline std::string toString(const T& v)
{
    std::ostringstream o;
    o << v;
    return o.str();
}

inline std::string toString(const std::string& s)
{
    return s;
}

template<class T> inline void fromString(const std::string& s, T& v)
{
    std::istringstream i(s);
    i >> v;
}

inline void fromString(const std::string& s, std::string& v)
{
    v = s;
}

template<class T>
class flat_data_trait<T, TYPE_SINGLE>
{
public:
    enum { is_struct = 0 };
    enum { is_vector = 0 };
    enum { is_single = 1 };
    typedef T data_type;
    typedef T value_type;
    static int size() { return 1; }
    static int size(const data_type&) { return size(); }
    static const char* header(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static const value_type* get(const data_type& d, int /*i*/ = 0) { return &d; }
    static void set(const value_type& v, data_type& d, int /*i*/ = 0) { d = v; }
    static void setS(const std::string& v, data_type& d, int /*i*/ = 0)
    {
        fromString(v, d);
    }
};

template<class T, int N = struct_data_trait<T>::NVAR>
class flat_struct_data_trait
{
public:
    enum { is_struct = 1 };
    enum { is_vector = 0 };
    enum { is_single = 0 };
    enum { is_default = ((struct_data_trait<T>::NVAR > 1) ? 1 : 0) };
    typedef T data_type;
    typedef std::string value_type;
    typedef struct_data_trait<data_type> shelper;
    typedef struct_data_trait_var<data_type,N-1> vhelper;
    typedef typename vhelper::value_type vtype;
    typedef default_flat_data_trait<vtype> vtrait;
    typedef typename vtrait::value_type iotype;
    typedef flat_struct_data_trait<data_type,N-1> prev;
    static int size() { return prev::size() + vtrait::size(); }
    static int size(const data_type&) { return size(); }
    static const char* header(const data_type& d, int i = 0)
    {
        int s = prev::size();
        if (i < s)
            return prev::header(d, i);
        else
        {
            const char* h1 = vhelper::shortname();
            const char* h2 = vtrait::header(*vhelper::get(d), i-s);
            if (h2 == NULL) return h1;
            else if (h1 == NULL) return h2;
            else
            {
                static std::string t;
                t = h1;
                t += " ";
                t += h2;
                return t.c_str();
            }
        }
    }
    static value_type* get(const data_type& d, int i = 0)
    {
        int s = prev::size();
        if (i < s)
            return prev::get(d, i);
        else
        {
            static std::string t;
            t = toString(*vtrait::get(*vhelper::get(d), i-s));
            return &t;
        }
    }
    static void setS(const std::string& v, data_type& d, int i = 0)
    {
        int s = prev::size();
        if (i < s)
            prev::setS(v, d, i);
        else
        {
            vtype var = *vhelper::get(d);
            vtrait::setS(v, var,i-s);
            vhelper::set(var, d);
        }
    }
    static void set(const value_type& v, data_type& d, int i = 0)
    {
        setS(v, d, i);
    }
};

template<class T>
class flat_struct_data_trait<T, 0>
{
public:
    enum { is_struct = 1 };
    enum { is_vector = 0 };
    enum { is_single = 0 };
    enum { is_default = ((struct_data_trait<T>::NVAR > 1) ? 1 : 0) };
    typedef T data_type;
    typedef std::string value_type;
    typedef struct_data_trait<data_type> shelper;
    static int size() { return 0; }
    static int size(const data_type&) { return size(); }
    static const char* header(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static value_type* get(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static void setS(const std::string& /*v*/, data_type& /*d*/, int /*i*/ = 0)
    {
    }
    static void set(const value_type& /*v*/, data_type& /*d*/, int /*i*/ = 0)
    {
    }
};

template<class T>
class flat_data_trait<T, TYPE_STRUCT> : public flat_struct_data_trait<T>
{
};

template<class T>
class flat_vector_data_trait
{
public:
    enum { is_struct = 0 };
    enum { is_vector = 1 };
    enum { is_single = 0 };
    enum { is_default = ((struct_data_trait<T>::NVAR == 1 && vector_data_trait<T>::NDIM >= 1) ? 1 : 0) };
    typedef T data_type;
    typedef vector_data_trait<data_type> vhelper;
    typedef typename vhelper::value_type vtype;
    typedef default_flat_data_trait<vtype> vtrait;
    typedef typename vtrait::value_type value_type;
    static int size() { return vhelper::SIZE * vtrait::size(); }
    static int size(const data_type&) { return size(); }
    static const char* header(const data_type& d, int i = 0)
    {
        int s = vtrait::size();
        int j = i / s;
        i = i % s;
        const char* h1 = vhelper::header(d, j);
        const char* h2 = vtrait::header(*vhelper::get(d, j), i);
        if (h2 == NULL) return h1;
        else if (h1 == NULL) return h2;
        else
        {
            static std::string t;
            t = h1;
            t += " ";
            t += h2;
            return t.c_str();
        }
    }
    static const value_type* get(const data_type& d, int i = 0)
    {
        int s = vtrait::size();
        int j = i / s;
        i = i % s;
        return vtrait::get(*vhelper::get(d, j), i);
    }
    static void set(const value_type& v, data_type& d, int i = 0)
    {
        int s = vtrait::size();
        int j = i / s;
        i = i % s;
        vtype t = *vhelper::get(d, j);
        vtrait::set(v, t, i);
        vhelper::set(t, d, j);
    }
    static void setS(const std::string& v, data_type& d, int i = 0)
    {
        int s = vtrait::size();
        int j = i / s;
        i = i % s;
        vtype t = *vhelper::get(d, j);
        vtrait::setS(v, t, i);
        vhelper::set(t, d, j);
    }
};

template<class T>
class flat_data_trait<T, TYPE_VECTOR> : public flat_vector_data_trait<T>
{
};

enum
{
    TABLE_NORMAL = 0,
    TABLE_HORIZONTAL = 1 << 0,
    TABLE_FIXEDSIZE = 1 << 1,
};

template<class T, int FLAGS = TABLE_NORMAL>
class table_data_widget_container
{
public:
    typedef T data_type;
    typedef vector_data_trait<data_type> rhelper;
    typedef typename rhelper::value_type row_type;
    //typedef vector_data_trait<row_type> vhelper;
    typedef default_flat_data_trait<row_type> vhelper;
    typedef typename vhelper::value_type value_type;
    typedef QVBoxLayout Layout;

    QSpinBox* wSize;
    QTableViewUpdater*  wTableView;
    QTableModelUpdater* wTableModel;
    QPushButtonUpdater* wDisplay;
    DataWidget * widget;
    Layout* container_layout;

    int rows;
    int cols;

    table_data_widget_container() : wSize(NULL), wTableView(NULL), wDisplay(NULL), widget(NULL), container_layout(NULL) {}


    bool createLayout( DataWidget* parent )
    {
        if( parent->layout() != NULL || container_layout != NULL) return false;
        container_layout = new Layout(parent);
        return true;
    }

    bool createLayout( QLayout* layout)
    {
        if ( container_layout != NULL ) return false;
        container_layout = new Layout(layout);
        return true;
    }


    void setCellText(int r, int c, const std::string& s)
    {
        if (FLAGS & TABLE_HORIZONTAL)
        {
            QStandardItem* item= wTableModel->item(c,r);
            if (!item) wTableModel->setItem(c, r, new QStandardItem(QString(s.c_str())));
            else       item->setText(QString(s.c_str()));
        }
        else
        {
            QStandardItem* item= wTableModel->item(r,c);
            if (!item) wTableModel->setItem(r, c, new QStandardItem(QString(s.c_str())));
            else       item->setText(QString(s.c_str()));
        }
    }
    std::string getCellText(int r, int c)
    {
        QStandardItem *item;
        QString text;
        if (FLAGS & TABLE_HORIZONTAL)
            item=wTableModel->item(c,r);
        else
            item=wTableModel->item(r,c);

        if (item)
        {
            text=item->text();

            if (!text.isNull())
                return std::string(text.ascii());
        }
        return std::string("");
    }

    template<class V>
    void setCell(int r, int c, const V& v)
    {
        std::ostringstream o;
        o << v;
        setCellText(r,c, o.str());
    }
    void setCell(int r, int c, const std::string& s)
    {
        setCellText(r, c, s);
    }
    template<class V>
    void getCell(int r, int c, V& v)
    {
        std::istringstream i(getCellText(r,c));
        i >> v;
    }
    void getCell(int r, int c, const std::string& s)
    {
        s = getCellText(r,c);
    }

    bool createWidgets(DataWidget* parent, const data_type& d, bool readOnly)
    {
        rows = 0;
        int dataRows = rhelper::size(d);

        wSize = new QSpinBox(0, INT_MAX, 1, parent);
		wSize->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);

        wDisplay = new QPushButtonUpdater( QString("Display the values"), parent);
		wDisplay->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);

        if (dataRows > 0)
            cols = vhelper::size(*rhelper::get(d,0));
        else
            cols = vhelper::size(row_type());


        if (FLAGS & TABLE_HORIZONTAL)
            wTableModel = new QTableModelUpdater(cols,dataRows,parent);
        else
            wTableModel = new QTableModelUpdater(dataRows,cols,parent);

        wTableView  = new QTableViewUpdater(parent);
        wTableView->setModel(wTableModel);

        widget=parent;

        wDisplay->setToggleButton(true);
        QWidget* parentWidget = parent;
        if(parentWidget)
            parentWidget = static_cast<QWidget*>(parentWidget->parent());
        bool propertyWidgetFlagOn = false;
        QDisplayDataWidget* displayDataWidget = dynamic_cast<QDisplayDataWidget*>(parentWidget);
        if(displayDataWidget)
            propertyWidgetFlagOn = displayDataWidget->flag().PROPERTY_WIDGET_FLAG;

		//if(isDisplayed() || propertyWidgetFlagOn)
		{
			processTableModifications(d);
			fillTable(d);
			rows = dataRows;
		}

        if(!propertyWidgetFlagOn)
            wDisplay->setOn(dataRows < MAX_NUM_ELEM && dataRows != 0 );
        wDisplay->setAutoDefault(false);

        wSize->setValue(dataRows);

        wTableView->setDisplayed(isDisplayed());

        if (readOnly)
        {
            wSize->setEnabled(false);
            wTableView->setEnabled(false);
        }
        else
        {
            if (!(FLAGS & TABLE_FIXEDSIZE))
            {
                parent->connect(wSize, SIGNAL( valueChanged(int) ), parent, SLOT( setWidgetDirty() ));
                //_widget->connect(wSize, SIGNAL( valueChanged(int) ), _widget, SLOT(updateDataValue()) );


                if( FLAGS & TABLE_HORIZONTAL)
                    parent->connect(wSize, SIGNAL( valueChanged(int) ), wTableModel, SLOT(resizeTableH(int) ) );
                else
                    parent->connect(wSize, SIGNAL( valueChanged(int) ), wTableModel, SLOT(resizeTableV(int) ) );
            }
            else
            {
                wSize->setEnabled(false);
            }
            parent->connect(wTableView, SIGNAL( activated(QModelIndex) ) , parent, SLOT(setWidgetDirty()) );
            parent->connect(wTableView, SIGNAL( pressed(QModelIndex) ) , parent, SLOT(setWidgetDirty()) );
            parent->connect(wTableView, SIGNAL( clicked(QModelIndex) ) , parent, SLOT(setWidgetDirty()) );
        }
        parent->connect(wDisplay, SIGNAL( toggled(bool) ), wTableView,   SLOT(setDisplayed(bool)));
        parent->connect(wDisplay, SIGNAL( toggled(bool) ), parent,   SLOT(setWidgetDirty()));
        parent->connect(wDisplay, SIGNAL( toggled(bool) ), wDisplay, SLOT(setDisplayed(bool)));
        parent->connect(wDisplay, SIGNAL( toggled(bool) ), parent, SLOT( updateWidgetValue() ));

        wDisplay->toggle();
        if(propertyWidgetFlagOn)
            wDisplay->setOn(false);
        else
            wDisplay->toggle();

        return true;
    }


    void setReadOnly(bool readOnly)
    {
        wSize->setEnabled(!readOnly);
        wTableView->setEnabled(!readOnly);
    }

    bool isDisplayed()
    {
        return (wDisplay->isOn());
    }


    void readFromData(const data_type& d)
    {
        if (isDisplayed())
        {
            processTableModifications(d);
            fillTable(d);
        }
    }
    void writeToData(data_type& d)
    {
        rows = wSize->value();
        if (!(FLAGS & TABLE_FIXEDSIZE))
        {
            int oldrows = rhelper::size(d);
            if( rows != oldrows)
            {
                rhelper::resize(rows,d);
            }
            int newrows = rhelper::size(d);
            if( rows != newrows)
            {
                wSize->setValue(newrows);
                rows = newrows;
                if (FLAGS & TABLE_HORIZONTAL)
                    wTableModel->setColumnCount(newrows);
                else
                    wTableModel->setRowCount(newrows);
            }
        }

        processTableModifications(d);

        if (isDisplayed())
        {
            for (int y=0; y<rows; ++y)
            {
                row_type r = *rhelper::get(d,y);
                for (int x=0; x<cols; ++x)
                {
                    value_type v = *vhelper::get(r,x);
                    getCell(y, x, v);
                    vhelper::set(v,r,x);
                }
                rhelper::set(r,d,y);
            }
        }

    }

    void processTableModifications(const data_type &d)
    {
        int currentTableNumRows;
        if (FLAGS & TABLE_HORIZONTAL)
            currentTableNumRows=wTableModel->columnCount();
        else
            currentTableNumRows=wTableModel->rowCount();
        int rows = wSize->value();

        QStringList horizontalHeaders;
        QStringList verticalHeaders;

        for (int x=0; x<cols; ++x)
        {
            const char* h = (rows > 0) ? vhelper::header(*rhelper::get(d,0),x) : vhelper::header(row_type(),x);

            if (h && *h)
                horizontalHeaders << h;
            else
                horizontalHeaders << QString::number(x);
        }

        if (FLAGS & TABLE_HORIZONTAL)
            wTableModel->setVerticalHeaderLabels(horizontalHeaders);
        else
            wTableModel->setHorizontalHeaderLabels(horizontalHeaders);


        for (int y=0; y<currentTableNumRows; ++y)
        {
            const char* h = rhelper::header(d,y);
            if (h && *h)
                verticalHeaders << h;
            else
                verticalHeaders << QString::number(y);
        }


        if (FLAGS & TABLE_HORIZONTAL)
            wTableModel->setHorizontalHeaderLabels(verticalHeaders);
        else
            wTableModel->setVerticalHeaderLabels(verticalHeaders);

        if (rows == currentTableNumRows)
        {
            return;
        }


        if (FLAGS & TABLE_HORIZONTAL)
            wTableModel->setColumnCount(rows);
        else
            wTableModel->setRowCount(rows);

        for (int y=currentTableNumRows; y<rows; ++y)
        {
            const char* h = rhelper::header(d,y);
            if (h && *h)
                verticalHeaders << h;
            else
                verticalHeaders << QString::number(y);
        }


        if (FLAGS & TABLE_HORIZONTAL)
            wTableModel->setHorizontalHeaderLabels(verticalHeaders);
        else
            wTableModel->setVerticalHeaderLabels(verticalHeaders);

    }

    void fillTable(const data_type &d)
    {
        int currentNum;
        int dsize = (int)d.size();
        if (FLAGS & TABLE_HORIZONTAL)
            currentNum=wTableModel->columnCount() > dsize ? dsize : wTableModel->columnCount();
        else
            currentNum=wTableModel->rowCount() > dsize ? dsize : wTableModel->rowCount();

        for (int y=0; y<currentNum; ++y)
            for (int x=0; x<cols; ++x)
                setCell(y, x, *vhelper::get(*rhelper::get(d,y),x));

    }

    void insertWidgets()
    {
        assert(container_layout);
        container_layout->add(wSize);
        container_layout->add(wTableView);
        container_layout->add(wDisplay);
    }
};

////////////////////////////////////////////////////////////////
/// variable-sized vectors support
////////////////////////////////////////////////////////////////

template<class T>
class vector_data_trait < std::vector<T> >
{
public:
    typedef std::vector<T> data_type;
    typedef T value_type;
    enum { NDIM = 1 };
    static int size(const data_type& d) { return d.size(); }
    static const char* header(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static const value_type* get(const data_type& d, int i = 0)
    {
        return ((unsigned)i < (unsigned)size(d)) ? &(d[i]) : NULL;
    }
    static void set( const value_type& v, data_type& d, int i = 0)
    {
        if ((unsigned)i < (unsigned)size(d))
            d[i] = v;
    }
    static void resize( int s, data_type& d)
    {
        d.resize(s);
    }
};




template<class T>
class vector_data_trait < sofa::helper::vector<T> > : public vector_data_trait< std::vector<T> >
{
};


// template<class T>
// class vector_data_trait < sofa::component::topology::PointData<T> > //: public vector_data_trait < sofa::helper::vector<T> >
// {
// public:
//     typedef sofa::component::topology::PointData<T> data_type;
//     typedef T value_type;
//     enum { NDIM = 1 };
//
//     static int size(const data_type& d) { return d.getValue().size(); }
//     static const char* header(const data_type& /*d*/, int /*i*/ = 0)
//     {
// 	return NULL;
//     }
//     static const value_type* get(const data_type& d, int i = 0)
//     {
// 	return ((unsigned)i < (unsigned)size(d.getValue())) ? &(d.getValue()[i]) : NULL;
//     }
//     static void set( const value_type& v, data_type& d, int i = 0)
//     {
// 	if ((unsigned)i < (unsigned)size(d.getValue()))
// 	{
// 	    sofa::helper::vector<T>& d_data = *(d.beginEdit());
// 	    d_data[i] = v;
// 	    d.endEdit();
// 	}
//     }
//     static void resize( int s, data_type& d)
//     {
//         sofa::helper::vector<T>& d_data = *(d.beginEdit());
// 	d_data.resize(s);
// 	d.endEdit();
//     }
// };

////////////////////////////////////////////////////////////////
/// sofa::defaulttype::ExtVector support
////////////////////////////////////////////////////////////////

template<class T>
class vector_data_trait < sofa::defaulttype::ExtVector<T> >
{
public:
    typedef sofa::defaulttype::ExtVector<T> data_type;
    typedef T value_type;
    enum { NDIM = 1 };
    static int size(const data_type& d) { return d.size(); }
    static const char* header(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static const value_type* get(const data_type& d, int i = 0)
    {
        return ((unsigned)i < (unsigned)size(d)) ? &(d[i]) : NULL;
    }
    static void set( const value_type& v, data_type& d, int i = 0)
    {
        if ((unsigned)i < (unsigned)size(d))
            d[i] = v;
    }
    static void resize( int s, data_type& d)
    {
        d.resize(s);
    }
};

template<class T>
class vector_data_trait < sofa::defaulttype::ResizableExtVector<T> > : public vector_data_trait < sofa::defaulttype::ExtVector<T> >
{
public:
};

////////////////////////////////////////////////////////////////
/// PointSubset support
////////////////////////////////////////////////////////////////
#ifdef TODOTOPO
template<>
class vector_data_trait < sofa::component::topology::PointSubset >
{
public:
    typedef sofa::component::topology::PointSubset data_type;
    typedef unsigned int value_type;
    enum { NDIM = 1 };
    static int size(const data_type& d) { return d.size(); }
    static const char* header(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static const value_type* get(const data_type& d, int i = 0)
    {
        return ((unsigned)i < (unsigned)size(d)) ? &(d[i]) : NULL;
    }
    static void set( const value_type& v, data_type& d, int i = 0)
    {
        if ((unsigned)i < (unsigned)size(d))
            d[i] = v;
    }
    static void resize( int s, data_type& d)
    {
        d.resize(s);
    }
};
#endif
////////////////////////////////////////////////////////////////
/// std::map from strings support
////////////////////////////////////////////////////////////////

template<class T>
class vector_data_trait < std::map<std::string, T> >
{
public:
    typedef std::map<std::string, T> data_type;
    typedef T value_type;
    enum { NDIM = 1 };
    static int size(const data_type& d) { return d.size(); }
    static const char* header(const data_type& d, int i = 0)
    {
        typename data_type::const_iterator it = d.begin();
        while (i > 0 && it != d.end())
        {
            ++it;
            --i;
        }
        if (i == 0) return it->first.c_str();
        else return NULL;
    }
    static const value_type* get(const data_type& d, int i = 0)
    {
        typename data_type::const_iterator it = d.begin();
        while (i > 0 && it != d.end())
        {
            ++it;
            --i;
        }
        if (i == 0) return &(it->second);
        else return NULL;
    }
    static void set( const value_type& v, data_type& d, int i = 0)
    {
        typename data_type::iterator it = d.begin();
        while (i > 0 && it != d.end())
        {
            ++it;
            --i;
        }
        if (i == 0) it->second = v;
    }
    static void resize( int /*s*/, data_type& /*d*/)
    {
        //d.resize(s);
    }
};



////////////////////////////////////////////////////////////////
/// deques support
////////////////////////////////////////////////////////////////

template<class T>
class vector_data_trait < std::deque<T> >
{
public:
    typedef std::deque<T> data_type;
    typedef T value_type;
    enum { NDIM = 1 };
    static int size(const data_type& d) { return d.size(); }
    static const char* header(const data_type& /*d*/, int /*i*/ = 0)
    {
        return NULL;
    }
    static const value_type* get(const data_type& d, int i = 0)
    {
        return ((unsigned)i < (unsigned)size(d)) ? &(d[i]) : NULL;
    }
    static void set( const value_type& v, data_type& d, int i = 0)
    {
        if ((unsigned)i < (unsigned)size(d))
            d[i] = v;
    }
    static void resize(int s, data_type& d)
    {
        d.resize(s);
    }
};




template<class T>
class vector_data_trait < sofa::helper::deque<T> > : public vector_data_trait< std::deque<T> >
{
};

} // namespace qt

} // namespace gui

} // namespace sofa


#endif
