/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
 *                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_SIMULATION_ANIMATION_LOOP_PARALLEL_SCHEDULER_H
#define SOFA_SIMULATION_ANIMATION_LOOP_PARALLEL_SCHEDULER_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/BaseAnimationLoop.h>
#include <sofa/core/ExecParams.h>
#include <sofa/simulation/common/common.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Visitor.h>
#include <sofa/helper/AdvancedTimer.h>

using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;



namespace sofa
{

namespace simulation
{


	class TaskScheduler;


class AnimationLoopParallelScheduler : public sofa::core::behavior::BaseAnimationLoop
{
public:

	typedef sofa::core::behavior::BaseAnimationLoop Inherit;
	SOFA_CLASS(AnimationLoopParallelScheduler,sofa::core::behavior::BaseAnimationLoop);

	Data<int> threadNumber;


protected:
	AnimationLoopParallelScheduler(simulation::Node* gnode = NULL);

	virtual ~AnimationLoopParallelScheduler();

public:
    virtual void init();

	/// Initialization method called at graph creation and modification, during bottom-up traversal.
	virtual void bwdInit();

	/// Update method called when variables used in precomputation are modified.
    virtual void reinit();

	virtual void cleanup();

	virtual void step(const core::ExecParams* params, double dt);


	/// Construction method called by ObjectFactory.
	template<class T>
    static typename T::SPtr create(T*, BaseContext* context, BaseObjectDescription* arg)
	{
		simulation::Node* gnode = dynamic_cast<simulation::Node*>(context);
		typename T::SPtr obj = sofa::core::objectmodel::New<T>(gnode);
		if (context) context->addObject(obj);
		if (arg) obj->parse(arg);
        return obj;
	}

private :

	// thread storage initialization
	void initThreadLocalData();

	int mNbThread;

	simulation::Node* gnode;
	


	//boost::shared_ptr<TaskScheduler> mScheduler;
};

} // namespace simulation

} // namespace sofa

#endif  /* SOFA_SIMULATION_ANIMATION_LOOP_PARALLEL_SCHEDULER_H */
