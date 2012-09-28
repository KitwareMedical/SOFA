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

#ifndef SOFA_COMPONENT_MAPPING_ARTICULATEDSYSTEMMAPPING_INL
#define SOFA_COMPONENT_MAPPING_ARTICULATEDSYSTEMMAPPING_INL

#include <sofa/component/mapping/ArticulatedSystemMapping.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/Multi2Mapping.inl>

#include <sofa/simulation/common/Simulation.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/gl/template.h>

#include <sofa/simulation/common/Node.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using sofa::simulation::Node;

using container::ArticulatedHierarchyContainer;
using container::ArticulationCenter;
using container::Articulation;

template <class TIn, class TInRoot, class TOut>
ArticulatedSystemMapping<TIn, TInRoot, TOut>::ArticulatedSystemMapping ()
    : ahc(NULL)
    , m_fromModel(NULL), m_toModel(NULL), m_fromRootModel(NULL)
{

}

template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::init()
{

    if(this->getFromModels1().empty())
    {
        serr << "Error while iniatilizing ; input Model not found" << sendl;
        return;
    }

    if(this->getToModels().empty())
    {
        serr << "Error while iniatilizing ; output Model not found" << sendl;
        return;
    }

    m_fromModel = this->getFromModels1()[0];
    m_toModel = this->getToModels()[0];

    Node* context = dynamic_cast<Node*>(m_fromModel->getContext());
    context->getNodeObject(ahc);
    articulationCenters = ahc->getArticulationCenters();

    const InVecCoord& xfrom = *m_fromModel->getX();

    ArticulationPos.clear();
    ArticulationAxis.clear();
    ArticulationPos.resize(xfrom.size());
    ArticulationAxis.resize(xfrom.size());

    //Root
    if(!this->getFromModels2().empty())
    {
        m_fromRootModel = this->getFromModels2()[0];
        sout << "Root Model found : Name = " << m_fromRootModel->getName() << sendl;
    }

    CoordinateBuf.clear();
    CoordinateBuf.resize(xfrom.size());
    for (unsigned int c=0; c<xfrom.size(); c++)
    {
        CoordinateBuf[c].x() = 0.0;
    }

    vector< ArticulationCenter* >::const_iterator ac = articulationCenters.begin();
    vector< ArticulationCenter* >::const_iterator acEnd = articulationCenters.end();

    for (; ac != acEnd; ac++)
    {
        (*ac)->OrientationArticulationCenter.clear();
        (*ac)->DisplacementArticulationCenter.clear();
        (*ac)->Disp_Rotation.clear();

        // sout << "(*ac)->OrientationArticulationCenter : " << (*ac)->OrientationArticulationCenter << sendl;
        // todo : warning if a (*a)->articulationIndex.getValue() exceed xfrom size !
    }

    helper::WriteAccessor<Data<OutVecCoord> > xtoData = *m_toModel->write(core::VecCoordId::position());
    apply(xtoData.wref(),
            xfrom,
            m_fromRootModel == NULL ? NULL : &m_fromRootModel->read(core::ConstVecCoordId::position())->getValue());
    Inherit::init();
    /*
    OutVecDeriv& vto = *m_toModel->getV();
    InVecDeriv& vfrom = *m_fromModel->getV();
    applyJT(vfrom, vto);
    */

}



template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::reset()
{
    init();
}



template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::apply( typename Out::VecCoord& out, const typename In::VecCoord& in, const typename InRoot::VecCoord* inroot  )
{
//    std::cout << " --> ArticulatedSystemMapping<TIn, TOut>::apply called with in: " << in << "  -- inroot" << (*inroot) << std::endl;

    const Data< OutVecCoord > &xtoData = *m_toModel->read(core::VecCoordId::position());
    out.resize(xtoData.getValue().size());

    // Copy the root position if a rigid root model is present
    if (m_fromRootModel && inroot)
    {
        out[0] = (*inroot)[m_fromRootModel->getSize()-1];
    }

    vector< ArticulationCenter* >::const_iterator ac = articulationCenters.begin();
    vector< ArticulationCenter* >::const_iterator acEnd = articulationCenters.end();

    for (; ac != acEnd; ac++)
    {
        int parent = (*ac)->parentIndex.getValue();
        int child = (*ac)->childIndex.getValue();

        // Before computing the child position, it is placed with the same orientation than its parent
        // and at the position compatible with the definition of the articulation center
        // (see initTranslateChild function for details...)
        Quat quat_child_buf = out[child].getOrientation();

        // The position of the articulation center can be deduced using the 6D position of the parent:
        // only useful for visualisation of the mapping => NO ! Used in applyJ and applyJT
        (*ac)->globalPosition.setValue(out[parent].getCenter() +
                out[parent].getOrientation().rotate((*ac)->posOnParent.getValue()));

        vector< Articulation* > articulations = (*ac)->getArticulations();
        vector< Articulation* >::const_iterator a = articulations.begin();
        vector< Articulation* >::const_iterator aEnd = articulations.end();

        int process = (*ac)->articulationProcess.getValue();

        switch(process)
        {
        case 0: // 0-(default) articulation are treated one by one, the axis of the second articulation is updated by the potential rotation of the first articulation
            //			   potential problems could arise when rotation exceed 90? (known problem of euler angles)
        {
            // the position of the child is reset to its rest position (based on the postion of the articulation center)
            out[child].getOrientation() = out[parent].getOrientation();
            out[child].getCenter() = out[parent].getCenter() + (*ac)->initTranslateChild(out[parent].getOrientation());

            Vec<3,OutReal> APos;
            APos = (*ac)->globalPosition.getValue();
            for (; a != aEnd; a++)
            {
                Vec<3,Real> axis = (*a)->axis.getValue();
                axis.normalize();
                (*a)->axis.setValue(axis);

                int ind = (*a)->articulationIndex.getValue();
                InCoord value = in[ind];
                axis = out[child].getOrientation().rotate((*a)->axis.getValue());
                ArticulationAxis[ind] = axis;

                if ((*a)->rotation.getValue())
                {
                    Quat dq;
                    dq.axisToQuat(axis, value.x());
                    out[child].getCenter() += (*ac)->translateChild(dq, out[child].getOrientation());
                    out[child].getOrientation() += dq;

                }
                if ((*a)->translation.getValue())
                {
                    out[child].getCenter() += axis*value.x();
                    APos += axis*value.x();
                }

                ArticulationPos[ind]= APos;
            }
            break;
        }
        case 1: // the axis of the articulations are linked to the parent - rotations are treated by successive increases -
        {
            //sout<<"Case 1"<<sendl;
            // no reset of the position of the child its position is corrected at the end to respect the articulation center.

            for (; a != aEnd; a++)
            {
                int ind = (*a)->articulationIndex.getValue();
                InCoord value = in[ind];
                InCoord prev_value = CoordinateBuf[ind];
                Vec<3,Real> axis = out[parent].getOrientation().rotate((*a)->axis.getValue());
                ArticulationAxis[ind]=axis;

                // the increment of rotation and translation are stored in dq and disp
                if ((*a)->rotation.getValue() )
                {
                    Quat r;
                    r.axisToQuat(axis, value.x() - prev_value.x());
                    // add the contribution into the quaternion that provides the actual orientation of the articulation center
                    (*ac)->OrientationArticulationCenter+=r;

                }
                if ((*a)->translation.getValue())
                {
                    (*ac)->DisplacementArticulationCenter+=axis*(value.x() - prev_value.x());
                }

            }

            //// in case 1: the rotation of the axis of the articulation follows the parent -> translation are treated "before":


            // step 1: compute the new position of the articulation center and the articulation pos
            //         rq: the articulation center folows the translations
            (*ac)->globalPosition.setValue(out[parent].getCenter() + out[parent].getOrientation().rotate((*ac)->posOnParent.getValue()) + (*ac)->DisplacementArticulationCenter);
            vector< Articulation* >::const_iterator a = articulations.begin();

            for (; a != aEnd; a++)
            {
                Vec<3,OutReal> APos;
                APos = (*ac)->globalPosition.getValue();
                ArticulationPos[(*a)->articulationIndex.getValue()]=APos;
            }

            // step 2: compute the position of the child
            out[child].getOrientation() = out[parent].getOrientation() + (*ac)->OrientationArticulationCenter;
            out[child].getCenter() =  (*ac)->globalPosition.getValue() - out[child].getOrientation().rotate( (*ac)->posOnChild.getValue() );

            break;

        }
        case 2: // the axis of the articulations are linked to the child (previous pos) - rotations are treated by successive increases -
        {
            //sout<<"Case 2"<<sendl;
            // no reset of the position of the child its position is corrected at the end to respect the articulation center.
            //Quat dq(0,0,0,1);
            Vec<3,Real> disp(0,0,0);

            for (; a != aEnd; a++)
            {
                int ind = (*a)->articulationIndex.getValue();
                InCoord value = in[ind];
                InCoord prev_value = CoordinateBuf[ind];
                Vec<3,Real> axis = quat_child_buf.rotate((*a)->axis.getValue());
                ArticulationAxis[ind]=axis;


                // the increment of rotation and translation are stored in dq and disp
                if ((*a)->rotation.getValue() )
                {
                    Quat r;
                    r.axisToQuat(axis, value.x() - prev_value.x());
                    // add the contribution into the quaternion that provides the actual orientation of the articulation center
                    (*ac)->OrientationArticulationCenter+=r;
                }
                if ((*a)->translation.getValue())
                {
                    disp += axis*(value.x()) ;

                }

                //// in case 2: the rotation of the axis of the articulation follows the child -> translation are treated "after"
                //// ArticulationPos do not move
                Vec<3,OutReal> APos;
                APos = (*ac)->globalPosition.getValue();
                ArticulationPos[(*a)->articulationIndex.getValue()]=APos;

            }
            (*ac)->DisplacementArticulationCenter=disp;

            out[child].getOrientation() = out[parent].getOrientation() + (*ac)->OrientationArticulationCenter;
            out[child].getCenter() =  (*ac)->globalPosition.getValue() - out[child].getOrientation().rotate((*ac)->posOnChild.getValue());
            out[child].getCenter() += (*ac)->DisplacementArticulationCenter;

            break;

        }
        }
    }

    //////////////////// buf the actual position of the articulations ////////////////////

    CoordinateBuf.clear();
    CoordinateBuf.resize(in.size());
    for (unsigned int c=0; c<in.size(); c++)
    {
        CoordinateBuf[c].x() = in[c].x();
    }

    //if( this->f_printLog.getValue())
    //{
    //	serr<<"ArticulatedSystemMapping::propageX processed :"<<sendl;
    //	if (m_fromRootModel!=NULL)
    //		serr<<"input root: "<<*m_fromRootModel->getX();
    //	serr<<"  - input: "<<*m_fromModel->getX()<<"  output : "<<*m_toModel->getX()<<sendl;
    //}

    //if( this->f_printLog.getValue())
    //{
    //	serr<<"ArticulatedSystemMapping::propageXfree processed"<<sendl;
    //	if (rootModel!=NULL)
    //		serr<<"input root: "<<*rootModel->getXfree();
    //	serr<<"  - input: "<<*m_fromModel->getXfree()<<"  output : "<<*m_toModel->getXfree()<<sendl;
    //}

//	  std::cout << " <-- ArticulatedSystemMapping<TIn, TOut>::apply called with in: " << in << "  -- inroot" << (*inroot) << std::endl;
}

template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in, const typename InRoot::VecDeriv* inroot )
{
    Data<OutVecCoord>* xtoData = m_toModel->write(core::VecCoordId::position());
    //const Data<InVecCoord>* xfromData = m_fromModel->read(core::ConstVecCoordId::position());

    /*apply(*xtoData->beginEdit(), xfromData->getValue(),
          (m_fromRootModel==NULL ? NULL : &m_fromRootModel->read(core::ConstVecCoordId::position())->getValue()));

    xtoData->endEdit();*/

    const OutVecCoord& xto = xtoData->getValue();

    /*std::cout << "--> applyJ : \n";
    std::cout << "xto=" << xto << std::endl;
    std::cout << "xfrom=" << xfromData->getValue() << std::endl;*/

    //sout<<" \n ApplyJ ";

    out.clear();
    out.resize(xto.size());

    // Copy the root position if a rigid root model is present
    if (m_fromRootModel && inroot)
    {
        // sout << "Root Model Name = " << rootModel->getName() << sendl;
        out[0] = (*inroot)[m_fromRootModel->getSize()-1];
    }
    else
        out[0] = OutDeriv();

    vector< ArticulationCenter* >::const_iterator ac = articulationCenters.begin();
    vector< ArticulationCenter* >::const_iterator acEnd = articulationCenters.end();

    int i = 0;

    for (; ac != acEnd; ac++)
    {
        int parent = (*ac)->parentIndex.getValue();
        int child = (*ac)->childIndex.getValue();

        getVOrientation(out[child]) += getVOrientation(out[parent]);
        Vec<3,OutReal> P = xto[parent].getCenter();
        Vec<3,OutReal> C = xto[child].getCenter();
        getVCenter(out[child]) = getVCenter(out[parent]) + cross(P-C, getVOrientation(out[parent]));
        //sout<<"P:"<< P  <<"- C: "<< C;

        vector< Articulation* > articulations = (*ac)->getArticulations();
        vector< Articulation* >::const_iterator a = articulations.begin();
        vector< Articulation* >::const_iterator aEnd = articulations.end();

        for (; a != aEnd; a++)
        {
            int ind = (*a)->articulationIndex.getValue();
            InCoord value = in[ind];
            Vec<3,OutReal> axis = ArticulationAxis[ind];
            Vec<3,OutReal> A = ArticulationPos[ind];


            if ((*a)->rotation.getValue())
            {
                getVCenter(out[child]) += cross(A-C, axis*value.x());
                getVOrientation(out[child]) += axis*value.x();
            }
            if ((*a)->translation.getValue())
            {
                getVCenter(out[child]) += axis*value.x();
            }
            i++;

        }
    }
    //if( this->f_printLog.getValue())
    //{
    //	serr<<" propagateV processed"<<sendl;
    //	if (m_fromRootModel!=NULL)
    //		serr<<"V input root: "<<*m_fromRootModel->getV();
    //	serr<<"  - V input: "<<*m_fromModel->getV()<<"   V output : "<<*m_toModel->getV()<<sendl;
    //}

    //if( this->f_printLog.getValue())
    //{
    //	serr<<"ArticulatedSystemMapping::propagateDx processed"<<sendl;
    //	if (m_fromRootModel!=NULL)
    //		serr<<"input root: "<<*m_fromRootModel->getDx();
    //	serr<<"  - input: "<<*m_fromModel->getDx()<<"  output : "<<*m_toModel->getDx()<<sendl;
    //}

    /*std::cout << "<-- applyJ : \n";
    std::cout << "xto=" << xto << std::endl;
    std::cout << "xfrom=" << xfromData->getValue() << std::endl;*/
}



template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in, typename InRoot::VecDeriv* outroot )
{
    //sout<<"\n ApplyJt";
    const OutVecCoord& xto = m_toModel->read(core::VecCoordId::position())->getValue();
//	InVecCoord &xfrom= *m_fromModel->read(core::ConstVecCoordId::position());

    //apply(xto,xfrom);

    // debug
    //apply(core::VecCoordId::position(), core::ConstVecCoordId::position(), (rootModel==NULL ? NULL : rootModel->getX()));
    //serr<<" XTO = "<<xto<<"  - Xroot :"<<*rootModel->getX()<<sendl;

    OutVecDeriv fObjects6DBuf = in;
    InVecDeriv OutBuf = out;

    vector< ArticulationCenter* >::const_iterator ac = articulationCenters.end();
    vector< ArticulationCenter* >::const_iterator acBegin = articulationCenters.begin();

    int i=ArticulationAxis.size();
    while (ac != acBegin)
    {
        ac--;
        int parent = (*ac)->parentIndex.getValue();
        int child = (*ac)->childIndex.getValue();

        getVCenter(fObjects6DBuf[parent]) += getVCenter(fObjects6DBuf[child]);
        Vec<3,OutReal> P = xto[parent].getCenter();
        Vec<3,OutReal> C = xto[child].getCenter();
        getVOrientation(fObjects6DBuf[parent]) += getVOrientation(fObjects6DBuf[child]) + cross(C-P,  getVCenter(fObjects6DBuf[child]));

        vector< Articulation* > articulations = (*ac)->getArticulations();

        vector< Articulation* >::const_iterator a = articulations.end();
        vector< Articulation* >::const_iterator aBegin = articulations.begin();

        while (a != aBegin)
        {
            a--;
            i--;
            int ind = (*a)->articulationIndex.getValue();
            Vec<3,OutReal> axis = ArticulationAxis[ind];
            Vec<3,Real> A = ArticulationPos[ind] ;
            OutDeriv T;
            getVCenter(T) = getVCenter(fObjects6DBuf[child]);
            getVOrientation(T) = getVOrientation(fObjects6DBuf[child]) + cross(C-A, getVCenter(fObjects6DBuf[child]));

            if ((*a)->rotation.getValue())
            {
                out[ind].x() += (InReal)dot(axis, getVOrientation(T));
            }
            if ((*a)->translation.getValue())
            {
                out[ind].x() += (InReal)dot(axis, getVCenter(T));
            }
        }
    }

    if (outroot && m_fromRootModel)
    {
        (*outroot)[m_fromRootModel->getSize()-1] += fObjects6DBuf[0];
    }

    //if( this->f_printLog.getValue())
    //{
    //	serr<<"ArticulatedSystemMapping::accumulateForce processed"<<sendl;
    //	serr<<" input f : "<<*m_toModel->getF();
    //	if (m_fromRootModel!=NULL)
    //		serr<<"- output root: "<<*m_fromRootModel->getF();
    //	serr<<"  - output F: "<<*m_fromModel->getF()<<sendl;
    //}

    //if( this->f_printLog.getValue())
    //{
    //	serr<<"ArticulatedSystemMapping::accumulateDf processed"<<sendl;
    //	serr<<" input df : "<<*m_toModel->getF();
    //	if (m_fromRootModel!=NULL)
    //		serr<<"- output root: "<<*m_fromRootModel->getF();
    //	serr<<"  - output: "<<*m_fromModel->getF()<<sendl;
    //}
}


template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::applyJT( InMatrixDeriv& out, const OutMatrixDeriv& in, InRootMatrixDeriv* outRoot )
{
    const OutVecCoord& xto = *m_toModel->getX();

    //std::cout << "applyJT (constraints) : \n";
    //std::cout << "xto = " << xto << std::endl;
    //std::cout << "xfrom = " << *m_fromModel->getX() << std::endl;
    //std::cout << "xfromFree = " << m_fromModel->read(core::VecCoordId::freePosition())->getValue() << std::endl;

    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt != colItEnd)
        {
            typename InMatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            //Hack to get a RowIterator, withtout default constructor
            InRootMatrixDeriv temp;
            typename InRootMatrixDeriv::RowIterator rootRowIt = temp.end();
            typename InRootMatrixDeriv::RowIterator rootRowItEnd = temp.end();

            if(m_fromRootModel && outRoot)
            {
                rootRowIt = outRoot->end();
                rootRowItEnd = outRoot->end();
            }

            while (colIt != colItEnd)
            {
                int childIndex = colIt.index();
                const OutDeriv valueConst = colIt.val();

                Vec<3,OutReal> C = xto[childIndex].getCenter();
                vector< ArticulationCenter* > ACList = ahc->getAcendantList(childIndex);

                vector< ArticulationCenter* >::const_iterator ac = ACList.begin();
                vector< ArticulationCenter* >::const_iterator acEnd = ACList.end();

                for (; ac != acEnd; ac++)
                {
                    vector< Articulation* > articulations = (*ac)->getArticulations();

                    vector< Articulation* >::const_iterator a = articulations.begin();
                    vector< Articulation* >::const_iterator aEnd = articulations.end();

                    for (; a != aEnd; a++)
                    {
                        int ind = (*a)->articulationIndex.getValue();
                        InDeriv data;

                        Vec< 3, OutReal > axis = ArticulationAxis[ind]; // xto[parent].getOrientation().rotate((*a)->axis.getValue());
                        Vec< 3, Real > A = ArticulationPos[ind] ; // Vec<3,OutReal> posAc = (*ac)->globalPosition.getValue();

                        OutDeriv T;
                        getVCenter(T) = getVCenter(valueConst);
                        getVOrientation(T) = getVOrientation(valueConst) + cross(C - A, getVCenter(valueConst));

                        if ((*a)->rotation.getValue())
                        {
                            data = (InReal)dot(axis, getVOrientation(T));
                        }

                        if ((*a)->translation.getValue())
                        {
                            data = (InReal)dot(axis, getVCenter(T));
                        }

                        o.addCol(ind, data);
                    }
                }

                if(m_fromRootModel && outRoot)
                {
                    unsigned int indexT = m_fromRootModel->getSize() - 1; // On applique sur le dernier noeud
                    Vec<3,OutReal> posRoot = xto[indexT].getCenter();

                    OutDeriv T;
                    getVCenter(T) = getVCenter(valueConst);
                    getVOrientation(T) = getVOrientation(valueConst) + cross(C - posRoot, getVCenter(valueConst));

                    if (rootRowIt == rootRowItEnd)
                        rootRowIt = (*outRoot).writeLine(rowIt.index());

                    rootRowIt.addCol(indexT, T);
                }

                ++colIt;
            }
        }
    }
}

template <class TIn, class TInRoot, class TOut>
void ArticulatedSystemMapping<TIn, TInRoot, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings()) return;
    std::vector< Vector3 > points;
    std::vector< Vector3 > pointsLine;

    vector< ArticulationCenter* >::const_iterator ac = articulationCenters.begin();
    vector< ArticulationCenter* >::const_iterator acEnd = articulationCenters.end();
    unsigned int i=0;
    for (; ac != acEnd; ac++)
    {
//		int parent = (*ac)->parentIndex.getValue();
//		int child = (*ac)->childIndex.getValue();
        vector< Articulation* > articulations = (*ac)->getArticulations();
        vector< Articulation* >::const_iterator a = articulations.begin();
        vector< Articulation* >::const_iterator aEnd = articulations.end();
        for (; a != aEnd; a++)
        {

            // Articulation Pos and Axis are based on the configuration of the parent
            int ind= (*a)->articulationIndex.getValue();
            points.push_back(ArticulationPos[ind]);

            pointsLine.push_back(ArticulationPos[ind]);
            Vec<3,OutReal> Pos_axis = ArticulationPos[ind] + ArticulationAxis[ind];
            pointsLine.push_back(Pos_axis);

            i++;
        }
    }

    vparams->drawTool()->drawPoints(points, 10, Vec<4,float>(1,0.5,0.5,1));
    vparams->drawTool()->drawLines(pointsLine, 1, Vec<4,float>(0,0,1,1));

    //
    //OutVecCoord& xto = *m_toModel->getX();
    //glDisable (GL_LIGHTING);
    //glPointSize(2);
}

} // namespace mapping

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_MAPPING_ARTICULATEDSYSTEMMAPPING_INL
