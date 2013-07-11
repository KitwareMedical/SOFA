#include "PythonMacros.h"

#include <sofa/core/objectmodel/Base.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/objectmodel/Context.h>
//#include <sofa/simulation/tree/GNode.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/core/BaseState.h>
#include <sofa/core/behavior/BaseMechanicalState.h>
#include <sofa/core/loader/BaseLoader.h>
#include <sofa/core/loader/MeshLoader.h>
#include <sofa/core/topology/Topology.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/MeshTopology.h>
#include <sofa/component/topology/GridTopology.h>
#include <sofa/component/topology/RegularGridTopology.h>
#include <sofa/component/typedef/Sofa_typedef.h>
#include <sofa/component/typedef/Mapping_double.h>
#include <sofa/core/BaseMapping.h>


#include "PythonScriptController.h"

#include "Binding_Base.h"
#include "Binding_BaseObject.h"
#include "Binding_BaseState.h"
#include "Binding_BaseMechanicalState.h"
#include "Binding_MechanicalObject.h"
#include "Binding_BaseContext.h"
#include "Binding_Context.h"
#include "Binding_Node.h"
#include "Binding_BaseLoader.h"
#include "Binding_MeshLoader.h"
#include "Binding_Topology.h"
#include "Binding_BaseMeshTopology.h"
#include "Binding_MeshTopology.h"
#include "Binding_GridTopology.h"
#include "Binding_RegularGridTopology.h"
#include "Binding_PythonScriptController.h"
#include "Binding_BaseMapping.h"
//#include "Binding_Mapping.h"
//#include "Binding_RigidMapping.h"
//#include "Binding_MultiMapping.h"
#include "Binding_SubsetMultiMapping.h"

// crée un objet Python à partir d'un objet Cpp héritant de Base,
// retournant automatiquement le type Python de plus haut niveau possible
// en fonction du type de l'objet Cpp
// Ceci afin de permettre l'utilisation de fonctions des sous-classes de Base
PyObject* SP_BUILD_PYSPTR(Base* obj)
{
    if (dynamic_cast<sofa::simulation::Node*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(Node));
    if (dynamic_cast<Context*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(Context));
    if (dynamic_cast<BaseContext*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseContext));

    if (dynamic_cast<sofa::core::loader::MeshLoader*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(MeshLoader));
    if (dynamic_cast<sofa::core::loader::BaseLoader*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseLoader));

    if (dynamic_cast<sofa::component::topology::RegularGridTopology*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(RegularGridTopology));
    if (dynamic_cast<sofa::component::topology::GridTopology*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(GridTopology));
    if (dynamic_cast<sofa::component::topology::MeshTopology*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(MeshTopology));
    if (dynamic_cast<sofa::component::topology::BaseMeshTopology*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseMeshTopology));
    if (dynamic_cast<sofa::component::topology::Topology*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(Topology));

    if (dynamic_cast<MechanicalObject3*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(MechanicalObject));
    if (dynamic_cast<sofa::core::behavior::BaseMechanicalState*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseMechanicalState));
    if (dynamic_cast<sofa::core::BaseState*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseState));

    if (dynamic_cast<sofa::component::controller::PythonScriptController*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(PythonScriptController));

    if (dynamic_cast<SubsetMultiMapping3d_to_3d*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(SubsetMultiMapping3_to_3));
    if (dynamic_cast<sofa::core::BaseMapping*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseMapping));

    if (dynamic_cast<BaseObject*>(obj))
        return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(BaseObject));

    // par défaut...
    return BuildPySPtr<Base>(obj,&SP_SOFAPYTYPEOBJECT(Base));
}
