#ifndef SOFA_COMPONENT_COLLISION_COMPLIANTCONTACT_H
#define SOFA_COMPONENT_COLLISION_COMPLIANTCONTACT_H

#include "BaseContact.h"
#include "../constraint/UnilateralConstraint.h"

#include "../initCompliant.h"

#include "../compliance/UniformCompliance.h"

#include "../mapping/ContactMapping.h" 		// should be normal mapping



namespace sofa
{

namespace component
{

namespace collision
{

template <class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes = sofa::defaulttype::Vec3Types>
class CompliantContact : public BaseContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>
{

public:

    SOFA_CLASS(SOFA_TEMPLATE3(CompliantContact, TCollisionModel1, TCollisionModel2, ResponseDataTypes), SOFA_TEMPLATE3(BaseContact, TCollisionModel1, TCollisionModel2, ResponseDataTypes) );

    typedef BaseContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes> Inherit;
    typedef typename Inherit::node_type node_type;
    typedef typename Inherit::delta_type delta_type;
    typedef typename Inherit::CollisionModel1 CollisionModel1;
    typedef typename Inherit::CollisionModel2 CollisionModel2;
    typedef typename Inherit::Intersection Intersection;
    Data< SReal > damping_ratio;
    Data< SReal > compliance_value;
    Data< SReal > restitution_coef;


protected:

    CompliantContact()
        : damping_ratio( initData(&damping_ratio, 0.0, "damping", "contact damping (use for stabilization)") )
        , compliance_value( initData(&compliance_value, 0.0, "compliance", "contact compliance: use model contact stiffnesses when < 0, use given value otherwise"))
        , restitution_coef( initData(&restitution_coef, 0.0, "restitution", "global restitution coef") )
    {}

    CompliantContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
        : Inherit(model1, model2, intersectionMethod)
        , damping_ratio( initData(&damping_ratio, (SReal)0.0, "damping", "contact damping (use for stabilization)") )
        , compliance_value( initData(&compliance_value, (SReal)0.0, "compliance", "contact compliance: use model contact stiffnesses when < 0, use given value otherwise"))
        , restitution_coef( initData(&restitution_coef, (SReal)0.0, "restitution", "global restitution coef") )
    {}

    typename node_type::SPtr create_node()
    {
        const unsigned size = this->mappedContacts.size();

        delta_type delta = this->make_delta();

        typename node_type::SPtr contact_node = node_type::create( this->getName() + " contact frame" );

        delta.node->addChild( contact_node.get() );

        // 1d contact dofs
        typedef container::MechanicalObject<defaulttype::Vec1Types> contact_dofs_type;
        typename contact_dofs_type::SPtr contact_dofs = sofa::core::objectmodel::New<contact_dofs_type>();

        contact_dofs->resize( size );
        contact_dofs->setName( this->getName() + " contact dofs" );
        contact_node->addObject( contact_dofs.get() );

        // contact mapping
        typedef mapping::ContactMapping<ResponseDataTypes, defaulttype::Vec1Types> contact_map_type;
        typename contact_map_type::SPtr contact_map = core::objectmodel::New<contact_map_type>();

        contact_map->setModels( delta.dofs.get(), contact_dofs.get() );
        contact_node->addObject( contact_map.get() );

        this->copyNormals( contact_map->normal );
        this->copyPenetrations( contact_map->penetrations );

        contact_map->init();	 

        // compliance
        typedef forcefield::UniformCompliance<defaulttype::Vec1Types> compliance_type;
        compliance_type::SPtr compliance = sofa::core::objectmodel::New<compliance_type>( contact_dofs.get() );
        contact_node->addObject( compliance.get() );
        compliance->compliance.setValue( compliance_value.getValue() );
        compliance->init();


        // projector
        typedef linearsolver::UnilateralConstraint projector_type;
        projector_type::SPtr projector = sofa::core::objectmodel::New<projector_type>(  );
        contact_node->addObject( projector.get() );
        
        // constraint value
        this->addConstraintValue( contact_node.get(), contact_dofs.get(), damping_ratio.getValue(), restitution_coef.getValue() );

        return delta.node;
    }

    void update_node(typename node_type::SPtr ) { }



};

} // namespace collision
} // namespace component
} // namespace sofa

#endif  // SOFA_COMPONENT_COLLISION_COMPLIANTCONTACT_H
