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

#include <sofa/component/collision/AddFramePerformer.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/collision/MouseInteractor.h>
#include <sofa/component/mapping/SkinningMapping.inl>
#include <sofa/helper/Quater.h>

namespace sofa
{

namespace component
{

namespace collision
{
template <class DataTypes>
void AddFramePerformer<DataTypes>::start()
{
    BodyPicked picked=this->interactor->getBodyPicked();
    if (!picked.body && !picked.mstate) return;

    vector<FBMapping*> vFBMapping;
    typename DataTypes::Coord point;

    if (picked.body)
    {
        point = picked.point;
        sofa::core::objectmodel::BaseContext* context=  picked.body->getContext();
        context->get<FBMapping>( &vFBMapping, core::objectmodel::BaseContext::SearchRoot);
    }
    else
    {
        core::behavior::MechanicalState<DataTypes>* mstateCollision=dynamic_cast< core::behavior::MechanicalState<DataTypes>*  >(picked.mstate);
        if (!mstateCollision)
        {
            this->interactor->serr << "incompatible MState during Mouse Interaction " << this->interactor->sendl;
            return;
        }
        static_cast<simulation::Node*>(mstateCollision->getContext())->get<FBMapping>( &vFBMapping, core::objectmodel::BaseContext::SearchRoot);
        int index = picked.indexCollisionElement;
        point=(*(mstateCollision->getX()))[index];
    }

    for( typename vector<FBMapping *>::iterator it = vFBMapping.begin(); it != vFBMapping.end(); it++)
        (*it)->insertFrame( point);
}

template <class DataTypes>
void AddFramePerformer<DataTypes>::execute()
{
};

template <class DataTypes>
AddFramePerformer<DataTypes>::AddFramePerformer(BaseMouseInteractor *i):TInteractionPerformer<DataTypes>(i)
{
}


template <class DataTypes>
AddFramePerformer<DataTypes>::~AddFramePerformer()
{
    //Should remove the frames added
};


}
}
}
