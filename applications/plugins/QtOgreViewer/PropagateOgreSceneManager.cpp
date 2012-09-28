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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "PropagateOgreSceneManager.h"
#include "OgreSceneObject.h"
#include <sofa/simulation/common/Node.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa
{
namespace simulation
{
namespace ogre
{

PropagateOgreSceneManager::PropagateOgreSceneManager(const core::ExecParams* params, Ogre::SceneManager* sceneMgr)
    :Visitor(params),
     mSceneMgr(sceneMgr)
{

}

Visitor::Result PropagateOgreSceneManager::processNodeTopDown(simulation::Node* node)
{
    Visitor::for_each(this, node, node->object, &PropagateOgreSceneManager::processObject);
    return RESULT_CONTINUE;
}

void PropagateOgreSceneManager::processObject(simulation::Node* /*node*/ ,core::objectmodel::BaseObject* obj)
{
    core::ogre::OgreSceneObject* ogreObj;
    if( (ogreObj = dynamic_cast<core::ogre::OgreSceneObject*>(obj)) != NULL )
    {
        ogreObj->setSceneManager(mSceneMgr);
    }
}
}
}
}
