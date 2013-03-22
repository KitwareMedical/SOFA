/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_VISUALMODEL_COLORMAP_H
#define SOFA_COMPONENT_VISUALMODEL_COLORMAP_H

#ifndef SOFA_NO_OPENGL

#include <sofa/core/objectmodel/Data.h>
#include <sofa/component/component.h>
#include <sofa/core/visual/VisualModel.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/rmath.h>
#include <sofa/helper/gl/template.h>
#include <sofa/defaulttype/Vec.h>
#include <string>

namespace sofa
{

namespace component
{

namespace visualmodel
{

class SOFA_OPENGL_VISUAL_API ColorMap : public virtual sofa::core::visual::VisualModel
{
public:
    SOFA_CLASS(ColorMap, sofa::core::visual::VisualModel);

    typedef defaulttype::Vec3f Color3;  // Color tripplet
    typedef defaulttype::Vec4f Color;   // ... with alpha value
    typedef sofa::helper::vector<Color> VecColor;

protected:
    ColorMap();
    virtual ~ColorMap() {
        glDeleteTextures(1, &texture);
    }

public:
    template<class Real>
    class evaluator
    {
    public:
        evaluator(const ColorMap* map, Real vmin, Real vmax)
            : map(map), vmin(vmin), vmax(vmax), vscale((vmax == vmin) ? (Real)0 : (map->entries.size()-1)/(vmax-vmin)) {}

        Color operator()(Real r)
        {
            Real e = (r-vmin)*vscale;
            if (e<0) return map->entries.front();

            unsigned int i = (unsigned int)(e);
            if (i>=map->entries.size()-1) return map->entries.back();

            Color c1 = map->entries[i];
            Color c2 = map->entries[i+1];
            return c1+(c2-c1)*(e-i);
        }
    protected:
        const ColorMap* map;
        const Real vmin;
        const Real vmax;
        const Real vscale;
    };

    Data<unsigned int> f_paletteSize;
    Data<sofa::helper::OptionsGroup> f_colorScheme;

    VecColor entries;
    GLuint texture;
    double min, max;

    void initOld(const std::string &data);

    void init();
    void reinit();

    //void initVisual() { initTextures(); }
    //void clearVisual() { }
    //void initTextures() {}
    void drawVisual(const core::visual::VisualParams* vparams);
    //void drawTransparent(const VisualParams* /*vparams*/)
    //void updateVisual();

    void prepareLegend();


    unsigned int getNbColors() { return entries.size(); }
    Color getColor(unsigned int i) {
        if (i < entries.size()) return entries[i];
        return Color(0.0, 0.0, 0.0, 0.0);
    }

    static ColorMap* getDefault();

    template<class Real>
    evaluator<Real> getEvaluator(Real vmin, Real vmax)
    {
        min = (double)vmin;
        max = (double)vmax;
        if (!entries.empty()) {
            return evaluator<Real>(this, vmin, vmax);
        } else {
            return evaluator<Real>(getDefault(), vmin, vmax);
        }
    }

    Color3 hsv2rgb(const Color3 &hsv);

    inline friend std::ostream& operator << (std::ostream& out, const ColorMap& m )
    {
        if (m.getName().empty()) out << "\"\"";
        else out << m.getName();
        out << " ";
        out << m.entries;
        return out;
    }

    inline friend std::istream& operator >> (std::istream& in, ColorMap& m )
    {
        std::string name;
        in >> name;
        if (name == "\"\"") m.setName("");
        else m.setName(name);
        in >> m.entries;
        return in;
    }
};

} // namespace visualmodel

} // namespace component

} // namespace sofa

#endif /* SOFA_NO_OPENGL */

#endif
