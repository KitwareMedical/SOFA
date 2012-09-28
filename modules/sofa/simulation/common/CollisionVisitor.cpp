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
#include <sofa/simulation/common/CollisionVisitor.h>
#include <sofa/simulation/common/Node.h>

namespace sofa
{

namespace simulation
{


Visitor::Result CollisionVisitor::processNodeTopDown(simulation::Node* node)
{
    for_each(this, node, node->collisionPipeline, &CollisionVisitor::processCollisionPipeline);
    return RESULT_CONTINUE;
}

void CollisionVisitor::processCollisionPipeline(simulation::Node*
#ifdef SOFA_DUMP_VISITOR_INFO
        node
#endif
        , core::collision::Pipeline* obj)
{
    //std::cerr<<"CollisionVisitor::processCollisionPipeline"<<std::endl;
#ifdef SOFA_DUMP_VISITOR_INFO
    printComment("computeCollisionReset");
    ctime_t t0=begin(node, obj);
#endif
    obj->computeCollisionReset();
#ifdef SOFA_DUMP_VISITOR_INFO
    end(node, obj,t0);
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    printComment("computeCollisionDetection");
    t0=begin(node, obj);
#endif
    obj->computeCollisionDetection();
#ifdef SOFA_DUMP_VISITOR_INFO
    end(node, obj,t0);
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    printComment("computeCollisionResponse");
    t0=begin(node, obj);
#endif
    obj->computeCollisionResponse();
#ifdef SOFA_DUMP_VISITOR_INFO
    end(node, obj,t0);
#endif
}

void CollisionResetVisitor::processCollisionPipeline(simulation::Node*
#ifdef SOFA_DUMP_VISITOR_INFO
        node
#endif
        , core::collision::Pipeline* obj)
{
#ifdef SOFA_DUMP_VISITOR_INFO
    printComment("computeCollisionReset");
    ctime_t t0=begin(node, obj);
#endif
    obj->computeCollisionReset();
#ifdef SOFA_DUMP_VISITOR_INFO
    end(node, obj,t0);
#endif
}

void CollisionDetectionVisitor::processCollisionPipeline(simulation::Node*
#ifdef SOFA_DUMP_VISITOR_INFO
        node
#endif
        , core::collision::Pipeline* obj)
{
#ifdef SOFA_DUMP_VISITOR_INFO
    printComment("computeCollisionDetection");
    ctime_t t0=begin(node, obj);
#endif
    obj->computeCollisionDetection();
#ifdef SOFA_DUMP_VISITOR_INFO
    end(node, obj,t0);
#endif
}

void CollisionResponseVisitor::processCollisionPipeline(simulation::Node*
#ifdef SOFA_DUMP_VISITOR_INFO
        node
#endif
        , core::collision::Pipeline* obj)
{
#ifdef SOFA_DUMP_VISITOR_INFO
    printComment("computeCollisionResponse");
    ctime_t t0=begin(node, obj);
#endif
    obj->computeCollisionResponse();
#ifdef SOFA_DUMP_VISITOR_INFO
    end(node, obj,t0);
#endif
}

#ifdef SOFA_SMP

void ParallelCollisionVisitor::processCollisionPipeline(simulation::Node*, core::collision::ParallelPipeline* obj)
{
    obj->parallelComputeCollisions();
}
void ParallelCollisionVisitor::processCollisionPipeline(simulation::Node* node, core::collision::Pipeline* obj)
{
    core::collision::ParallelPipeline* p = dynamic_cast<core::collision::ParallelPipeline*>(obj);
    if (p)
        processCollisionPipeline(node, p);
}

#endif

} // namespace simulation

} // namespace sofa

