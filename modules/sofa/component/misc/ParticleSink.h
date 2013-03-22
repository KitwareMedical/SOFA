/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Interface: ParticleSink
//
// Description:
//
//
// Author: Jeremie Allard, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_MISC_PARTICLESINK_H
#define SOFA_COMPONENT_MISC_PARTICLESINK_H

#include <sofa/helper/system/config.h>
#include <sofa/helper/gl/template.h>
#include <sofa/core/behavior/ProjectiveConstraintSet.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/component/topology/TopologySubsetData.inl>
#include <sofa/component/topology/PointSetTopologyModifier.h>
#include <sofa/core/topology/TopologyChange.h>
#include <sofa/component/container/MechanicalObject.h>
#include <vector>
#include <iterator>
#include <iostream>
#include <ostream>
#include <algorithm>

namespace sofa
{

namespace component
{

namespace misc
{

using namespace sofa::component::topology;

template<class TDataTypes>
class ParticleSink : public core::behavior::ProjectiveConstraintSet<TDataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ParticleSink,TDataTypes), SOFA_TEMPLATE(core::behavior::ProjectiveConstraintSet,TDataTypes));

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowType MatrixDerivRowType;
    typedef helper::vector<Real> VecDensity;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalModel;
    typedef helper::vector<unsigned int> SetIndexArray;

    Data<Deriv> planeNormal;
    Data<Real> planeD0;
    Data<Real> planeD1;
    Data<defaulttype::Vec3f> color;
    Data<bool> showPlane;

    sofa::component::topology::PointSubsetData< SetIndexArray > fixed;
    //Data< SetIndexArray > fixed;
protected:
    ParticleSink()
        : planeNormal(initData(&planeNormal, "normal", "plane normal"))
        , planeD0(initData(&planeD0, (Real)0, "d0", "plane d coef at which particles acceleration is constrained to 0"))
        , planeD1(initData(&planeD1, (Real)0, "d1", "plane d coef at which particles are removed"))
        , color(initData(&color, defaulttype::Vec3f(0.0f,.5f,.2f), "color", "plane color"))
        , showPlane(initData(&showPlane, false, "showPlane", "enable/disable drawing of plane"))
        , fixed(initData(&fixed, "fixed", "indices of fixed particles"))
    {
        this->f_listening.setValue(true);
        Deriv n;
        DataTypes::set(n, 0, 1, 0);
        planeNormal.setValue(n);
    }

    virtual ~ParticleSink()
    {
    }
public:
    virtual void init()
    {
        this->core::behavior::ProjectiveConstraintSet<TDataTypes>::init();
        if (!this->mstate) return;

        sout << "ParticleSink: normal="<<planeNormal.getValue()<<" d0="<<planeD0.getValue()<<" d1="<<planeD1.getValue()<<sendl;

        sofa::core::topology::BaseMeshTopology* _topology;
        _topology = this->getContext()->getMeshTopology();

        // Initialize functions and parameters for topology data and handler
        fixed.createTopologicalEngine(_topology);
        fixed.registerTopologicalData();

    }

    virtual void animateBegin(double /*dt*/, double time)
    {
        //sout << "ParticleSink: animate begin time="<<time<<sendl;
        if (!this->mstate) return;
        const VecCoord& x = *this->mstate->getX();
        const VecDeriv& v = *this->mstate->getV();
        int n = x.size();
        helper::vector<unsigned int> remove;
        const bool log = this->f_printLog.getValue();
        for (int i=n-1; i>=0; --i) // always remove points in reverse order
        {
            Real d = x[i]*planeNormal.getValue()-planeD1.getValue();
            if (d<0)
            {
                if (log)
                    sout << "SINK particle "<<i<<" time "<<time<<" position "<<x[i]<<" velocity "<<v[i]<<sendl;
                remove.push_back(i);
            }
        }
        if (!remove.empty())
        {
            //sofa::core::topology::BaseMeshTopology* _topology;
            //_topology = this->getContext()->getMeshTopology();

            sofa::component::topology::PointSetTopologyModifier* pointMod;
            this->getContext()->get(pointMod);

            if (pointMod != NULL)
            {
                sout << "ParticleSink: remove "<<remove.size()<<" particles using PointSetTopologyModifier."<<sendl;
                pointMod->removePointsWarning(remove);
                pointMod->propagateTopologicalChanges();
                pointMod->removePointsProcess(remove);
            }
            else if(container::MechanicalObject<DataTypes>* object = dynamic_cast<container::MechanicalObject<DataTypes>*>(this->mstate.get()))
            {
                sout << "ParticleSink: remove "<<remove.size()<<" particles using MechanicalObject."<<sendl;
                // deleting the vertices
                for (unsigned int i = 0; i < remove.size(); ++i)
                {
                    --n;
                    object->replaceValue(n, remove[i] );
                }
                // resizing the state vectors
                this->mstate->resize(n);
            }
            else
            {
                sout << "ERROR(ParticleSink): no external object supporting removing points!"<<sendl;
            }
        }
    }


    template <class DataDeriv>
    void projectResponseT(DataDeriv& res) ///< project dx to constrained space
    {
        if (!this->mstate) return;
        if (fixed.getValue().empty()) return;

        const SetIndexArray& _fixed = fixed.getValue();
        // constraint the last value
        for (unsigned int s=0; s<_fixed.size(); s++)
            res[_fixed[s]] = Deriv();
    }
    virtual void projectResponse(VecDeriv& res) ///< project dx to constrained space
    {
        projectResponseT(res);
    }
    virtual void projectResponse(MatrixDerivRowType& res) ///< project dx to constrained space
    {
        projectResponseT(res);
    }

    virtual void projectVelocity(VecDeriv&) ///< project dx to constrained space (dx models a velocity)
    {
    }

    virtual void projectPosition(VecCoord& x) ///< project x to constrained space (x models a position)
    {
        if (!this->mstate) return;

        helper::WriteAccessor< Data< SetIndexArray > > _fixed = fixed;

        _fixed.clear();
        // constraint the last value
        for (unsigned int i=0; i<x.size(); i++)
        {
            Real d = x[i]*planeNormal.getValue()-planeD0.getValue();
            if (d<0)
            {
                _fixed.push_back(i);
            }
        }
    }

    virtual void animateEnd(double /*dt*/, double /*time*/)
    {

    }

    virtual void handleEvent(sofa::core::objectmodel::Event* event)
    {
        if (simulation::AnimateBeginEvent* ev = dynamic_cast<simulation::AnimateBeginEvent*>(event))
            animateBegin(ev->getDt(), this->getContext()->getTime());
        if (simulation::AnimateEndEvent* ev = dynamic_cast<simulation::AnimateEndEvent*>(event))
            animateEnd(ev->getDt(), this->getContext()->getTime());
    }

    virtual void draw()
    {
#ifndef SOFA_NO_OPENGL
        if (!showPlane.getValue()) return;
        defaulttype::Vec3d normal; normal = planeNormal.getValue();

        // find a first vector inside the plane
        defaulttype::Vec3d v1;
        if( 0.0 != normal[0] ) v1 = defaulttype::Vec3d(-normal[1]/normal[0], 1.0, 0.0);
        else if ( 0.0 != normal[1] ) v1 = defaulttype::Vec3d(1.0, -normal[0]/normal[1],0.0);
        else if ( 0.0 != normal[2] ) v1 = defaulttype::Vec3d(1.0, 0.0, -normal[0]/normal[2]);
        v1.normalize();
        // find a second vector inside the plane and orthogonal to the first
        defaulttype::Vec3d v2;
        v2 = v1.cross(normal);
        v2.normalize();
        const float size=1.0f;
        defaulttype::Vec3d center = normal*planeD0.getValue();
        defaulttype::Vec3d corners[4];
        corners[0] = center-v1*size-v2*size;
        corners[1] = center+v1*size-v2*size;
        corners[2] = center+v1*size+v2*size;
        corners[3] = center-v1*size+v2*size;

        // glEnable(GL_LIGHTING);
        glDisable(GL_LIGHTING);
        glEnable(GL_CULL_FACE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glCullFace(GL_FRONT);

        glColor3f(color.getValue()[0],color.getValue()[1],color.getValue()[2]);

        glBegin(GL_QUADS);
        helper::gl::glVertexT(corners[0]);
        helper::gl::glVertexT(corners[1]);
        helper::gl::glVertexT(corners[2]);
        helper::gl::glVertexT(corners[3]);
        glEnd();

        glDisable(GL_CULL_FACE);

        glColor4f(1,0,0,1);
#endif /* SOFA_NO_OPENGL */
    }
};

}

} // namespace component

} // namespace sofa

#endif

