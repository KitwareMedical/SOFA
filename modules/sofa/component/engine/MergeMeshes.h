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
#ifndef SOFA_COMPONENT_ENGINE_MERGEMESHS_H
#define SOFA_COMPONENT_ENGINE_MERGEMESHS_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/loader/MeshLoader.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * This class merge several meshes.
 */
template <class DataTypes>
class MergeMeshes : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(MergeMeshes,DataTypes),core::DataEngine);
    typedef typename DataTypes::VecCoord VecCoord;
    typedef helper::vector<unsigned int> VecIndex;

protected:
    MergeMeshes();

    ~MergeMeshes();
public:
    /// Parse the given description to assign values to this object's fields and potentially other parameters
    void parse ( sofa::core::objectmodel::BaseObjectDescription* arg );

    /// Assign the field values stored in the given map of name -> value pairs
    void parseFields ( const std::map<std::string,std::string*>& str );

    void init();

    void reinit();

    void update();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const MergeMeshes<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    Data<unsigned int> f_nbMeshes;

    helper::vector<Data<VecCoord>*> vf_positions;
    helper::vector< Data< helper::vector< helper::fixed_array <unsigned int,2> > >* > vf_edges;
    helper::vector< Data< helper::vector< helper::fixed_array <unsigned int,3> > >* > vf_triangles;
    helper::vector< Data< helper::vector< helper::fixed_array <unsigned int,4> > >* > vf_quads;
    helper::vector< Data< helper::vector< helper::vector <unsigned int> > >* > vf_polygons;
    helper::vector< Data< helper::vector< helper::fixed_array<unsigned int,4> > >* > vf_tetrahedra;
    helper::vector< Data< helper::vector< helper::fixed_array<unsigned int,8> > >* > vf_hexahedra;


    Data<unsigned> f_output_npoints;
    Data<VecCoord> f_output_positions;
    Data< helper::vector< helper::fixed_array <unsigned int,2> > > f_output_edges;
    Data< helper::vector< helper::fixed_array <unsigned int,3> > > f_output_triangles;
    Data< helper::vector< helper::fixed_array <unsigned int,4> > > f_output_quads;
    Data< helper::vector< helper::vector <unsigned int> > > f_output_polygons;
    Data< helper::vector< helper::fixed_array<unsigned int,4> > > f_output_tetrahedra;
    Data< helper::vector< helper::fixed_array<unsigned int,8> > > f_output_hexahedra;

protected:
    void createInputMeshesData(int nb = -1);
    template<class T>
    void createInputDataVector(unsigned int nb, helper::vector< Data<T>* >& vf, std::string name, std::string help)
    {
        vf.reserve(nb);
        for (unsigned int i=vf.size(); i<nb; ++i)
        {
            std::ostringstream oname, ohelp;
            oname << name << (i+1);
            ohelp << help << (i+1);
            std::string name_i = oname.str();
            std::string help_i = ohelp.str();
            Data<T>* d = new Data<T>(help_i.c_str(), true, false);
            d->setName(name_i);
            vf.push_back(d);
            this->addData(d);
            this->addInput(d);
        }
    }
    template<class T>
    void deleteInputDataVector(helper::vector< Data<T>* >& vf)
    {
        for (unsigned int i=0; i<vf.size(); ++i)
        {
            this->delInput(vf[i]);
            delete vf[i];
        }
        vf.clear();
    }
    template<class T>
    void mergeInputDataVector(unsigned int nb, Data<T>& outF, const helper::vector< Data<T>* >& inVF, const helper::vector< Data<VecCoord>* >& inVFPos)
    {
        unsigned int nelems = 0;
        for (unsigned int i=0; i<nb; ++i)
            nelems += inVF[i]->getValue().size();
        helper::WriteAccessor< Data<T> > out = outF;
        out.clear();
        out.reserve(nelems);
        unsigned int shift = 0;
        for (unsigned int i=0; i<nb; ++i)
        {
            helper::ReadAccessor< Data<T> > in = inVF[i];
            for (unsigned int j=0; j<in.size(); ++j)
            {
                typename T::value_type outT = in[j];
                for (unsigned int k=0; k<outT.size(); ++k)
                    outT[k] += shift;
                out.push_back(outT);
            }
            shift += inVFPos[i]->getValue().size();
        }
    }
    void mergeInputDataVector(unsigned int nb, Data<VecCoord>& outF, const helper::vector< Data<VecCoord>* >& inVF)
    {
        unsigned int nelems = 0;
        for (unsigned int i=0; i<nb; ++i)
            nelems += inVF[i]->getValue().size();
        helper::WriteAccessor< Data<VecCoord> > out = outF;
        out.clear();
        out.reserve(nelems);
        for (unsigned int i=0; i<nb; ++i)
        {
            helper::ReadAccessor< Data<VecCoord> > in = inVF[i];
            for (unsigned int j=0; j<in.size(); ++j)
                out.push_back(in[j]);
        }
    }
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_MERGEMESHES_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Vec1dTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Vec2dTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Vec3dTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Rigid2dTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Vec1fTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Vec2fTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Vec3fTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Rigid2fTypes>;
extern template class SOFA_ENGINE_API MergeMeshes<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
