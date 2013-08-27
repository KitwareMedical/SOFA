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
#ifndef SOFA_COMPONENT_FORCEFIELD_STANDARDTETRAHEDRALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_STANDARDTETRAHEDRALFEMFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif
#include <sofa/component/fem/material/HyperelasticMaterial.h>
#include <sofa/component/component.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/MatSym.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/TopologyData.h>
#include <string>
#include <map>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace forcefield
{
using namespace std;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;


//***************** Tetrahedron FEM code for several elastic models: StandardTetrahedralFEMForceField*******************************************************************
//********************************** Based on classical discretization : Fi=-Bi^T S V and Kij=Bi^T N Bj +Di^T S Dj **********************************************
//***************************************** where Bi is the strain displacement (6*3 matrix), S SPK tensor N=dS/dC, Di shape vector ************************************
//**************************** Code dependant on HyperelasticMatrialFEM and inherited classes *********************************************************************

/** Compute Finite Element forces based on tetrahedral elements.
*/
template<class DataTypes>
class StandardTetrahedralFEMForceField: public core::behavior::ForceField<DataTypes>
{
  public:
	  SOFA_CLASS(SOFA_TEMPLATE(StandardTetrahedralFEMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef Mat<3,3,Real> Matrix3;
	typedef Mat<6,6,Real> Matrix6;
	typedef Mat<6,3,Real> Matrix63;
	typedef MatSym<3,Real> MatrixSym;

	typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv; 
	typedef core::objectmodel::Data<VecCoord>    DataVecCoord; 

    typedef helper::vector<Real> SetParameterArray;
    typedef helper::vector<Coord> SetAnisotropyDirectionArray;
	

    typedef core::topology::BaseMeshTopology::index_type Index;
    typedef core::topology::BaseMeshTopology::Tetra Element;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecElement;



public :
	
	typename sofa::component::fem::MaterialParameters<DataTypes> globalParameters;

	
	

    /// data structure stored for each tetrahedron
	class TetrahedronRestInformation : public fem::StrainInformation<DataTypes>
    {
    public:

      /// rest volume
	  Real restVolume;
      /// current tetrahedron volume
      Real volScale;

      /// shape vector at the rest configuration
	  Coord shapeVector[4];

	  /// fiber direction in rest configuration
      Coord fiberDirection;
	 
	  /// derivatives of J
	  Coord dJ[4];
	  MatrixSym SPKTensorGeneral;
	  /// deformation gradient = gradPhi
	  Matrix3 deformationGradient;
	  /// right Cauchy-Green deformation tensor C (gradPhi^T gradPhi) 
	  Matrix63 matB[4];
	  Real strainEnergy;

      //Tetrahedron Points Indicies for CUDA
      float tetraIndices[4];
      //Tetrahedron Edges for CUDA
      float tetraEdges[6];

      /// Output stream
      inline friend ostream& operator<< ( ostream& os, const TetrahedronRestInformation& /*eri*/ ) {  return os;  }
      /// Input stream
      inline friend istream& operator>> ( istream& in, TetrahedronRestInformation& /*eri*/ ) { return in; }

      TetrahedronRestInformation() {}  
    };
    typedef typename VecCoord::template rebind<TetrahedronRestInformation>::other tetrahedronRestInfoVector;
    
	
   /// data structure stored for each edge
   class EdgeInformation
   {
   public:
	   /// store the stiffness edge matrix 
	   Matrix3 DfDx;
       float vertices[2];

	   /// Output stream
	   inline friend ostream& operator<< ( ostream& os, const EdgeInformation& /*eri*/ ) {  return os;  }
	   /// Input stream
	   inline friend istream& operator>> ( istream& in, EdgeInformation& /*eri*/ ) { return in; }

     EdgeInformation() {}
   };
   typedef typename VecCoord::template rebind<EdgeInformation>::other edgeInformationVector;

 protected :
   core::topology::BaseMeshTopology* _topology;
   VecCoord  _initialPoints;	/// the intial positions of the points
   bool updateMatrix;
   bool  _meshSaved ;
   Data<string> f_materialName; /// the name of the material
   Data<SetParameterArray> f_parameterSet;
   Data<SetAnisotropyDirectionArray> f_anisotropySet;
   Data<string> f_parameterFileName;

   
public:

	void setMaterialName(const string name) {
		f_materialName.setValue(name);
	}
	void setparameter(const vector<Real> param) {
		f_parameterSet.setValue(param);
	}
	void setdirection(const vector<Coord> direction) {
		f_anisotropySet.setValue(direction);
	}

protected:
   StandardTetrahedralFEMForceField();
   
   virtual   ~StandardTetrahedralFEMForceField();
public:

  //  virtual void parse(core::objectmodel::BaseObjectDescription* arg);

    virtual void init();
    //Used for CUDA implementation
    void initNeighbourhoodPoints();
    void initNeighbourhoodEdges();
    
    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v);
    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx);

    void draw(const core::visual::VisualParams* vparams);

	Mat<3,3,double> getPhi( int );

        class GHTetrahedronHandler : public TopologyDataHandler<Tetrahedron, tetrahedronRestInfoVector >
        {
        public:
          typedef typename StandardTetrahedralFEMForceField<DataTypes>::TetrahedronRestInformation TetrahedronRestInformation;

          GHTetrahedronHandler(StandardTetrahedralFEMForceField<DataTypes>* ff,
                                 TetrahedronData<tetrahedronRestInfoVector>* data )
            :TopologyDataHandler<Tetrahedron, tetrahedronRestInfoVector >(data)
            ,ff(ff)
          {
          }

          void applyCreateFunction(unsigned int, TetrahedronRestInformation &t, const Tetrahedron
                                   &, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);

         protected:
          StandardTetrahedralFEMForceField<DataTypes>* ff;

        };
	
  protected:
    /// the array that describes the complete material energy and its derivatives

	fem::HyperelasticMaterial<DataTypes> *myMaterial;

        TetrahedronData<tetrahedronRestInfoVector> tetrahedronInfo;
        //EdgeData<sofa::helper::vector< EdgeInformation> > edgeInfo;
        EdgeData<edgeInformationVector> edgeInfo;


        void testDerivatives();
        void saveMesh( const char *filename );
	
	VecCoord myposition;

        GHTetrahedronHandler* tetrahedronHandler;
    
};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif
#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_STANDARDTETRAHEDRALFEMFORCEFIELD_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_MISC_FEM_API StandardTetrahedralFEMForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_FEM_API StandardTetrahedralFEMForceField<Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_STANDARDTETRAHEDRALFEMFORCEFIELD_CPP)


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_STANDARDTETRAHEDRALFEMFORCEFIELD_H
