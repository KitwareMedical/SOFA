/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* version 0.1                                                                  *
* Copyright (C) 2009-2012, IGG Team, LSIIT, University of Strasbourg           *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#ifndef __IMPLICIT_HIERARCHICAL_MAP3__
#define __IMPLICIT_HIERARCHICAL_MAP3__

#include "Topology/map/embeddedMap3.h"

namespace CGoGN
{

namespace Algo
{

namespace Volume
{

namespace IHM
{


template<typename T, unsigned int ORBIT> class AttributeHandler_IHM ;

class ImplicitHierarchicalMap3 : public EmbeddedMap3
{
	template<typename T, unsigned int ORBIT> friend class AttributeHandler_IHM ;

public:

	FunctorType* vertexVertexFunctor ;
	FunctorType* edgeVertexFunctor ;
	FunctorType* faceVertexFunctor ;
	FunctorType* volumeVertexFunctor ;

	unsigned int m_curLevel ;
	unsigned int m_maxLevel ;
	unsigned int m_edgeIdCount ;
	unsigned int m_faceIdCount;

	DartAttribute<unsigned int> m_dartLevel ;
	DartAttribute<unsigned int> m_edgeId ;
	DartAttribute<unsigned int> m_faceId ;

	AttributeMultiVector<unsigned int>* m_nextLevelCell[NB_ORBITS] ;

public:
	ImplicitHierarchicalMap3() ;

	~ImplicitHierarchicalMap3() ;

	static const unsigned int DIMENSION = 3 ;

	//!
	/*!
	 *
	 */
	void update_topo_shortcuts();

	//!
	/*!
	 *
	 */
	void initImplicitProperties() ;

	/**
	 * clear the map
	 * @param remove attrib remove attribute (not only clear the content)
	 */
	void clear(bool removeAttrib);

	/*! @name Attributes Management
	 *  To handles Attributes for each level of an implicit 3-map
	 *************************************************************************/

	//@{
	//!
	/*!
	 *
	 */
	template <typename T, unsigned int ORBIT>
	AttributeHandler_IHM<T, ORBIT> addAttribute(const std::string& nameAttr) ;

	//!
	/*!
	 *
	 */
	template <typename T, unsigned int ORBIT>
	AttributeHandler_IHM<T, ORBIT> getAttribute(const std::string& nameAttr) ;
	//@}

	/*! @name Basic Topological Operators
	 *  Redefinition of the basic topological operators
	 *************************************************************************/

	//@{
	virtual Dart newDart() ;

	Dart phi1(Dart d) ;

	Dart phi_1(Dart d) ;

	Dart phi2(Dart d) ;

private:
	Dart phi2bis(Dart d) ;

public:
	Dart phi3(Dart d);

	Dart alpha0(Dart d);

	Dart alpha1(Dart d);

	Dart alpha2(Dart d);

	Dart alpha_2(Dart d);
	//@}

	/*! @name Topological Operators with Cells id management
	 *  Topological operations on Hierarchical Implicit 3-maps
	 *************************************************************************/

//	void deleteVolume(Dart d);

	bool isWellEmbedded();

	//@{

	//!
	/*!
	 *
	 */
	void swapEdges(Dart d, Dart e);

	//!
	/*!
	 *
	 */
	void saveRelationsAroundVertex(Dart d, std::vector<std::pair<Dart, Dart> >& vd);

	void unsewAroundVertex(std::vector<std::pair<Dart, Dart> >& vd);

	Dart quadranguleFace(Dart d);

	void deleteVertexSubdividedFace(Dart d);
	//@}

	void setVertexVertexFunctor(FunctorType* f) { vertexVertexFunctor = f ; }
	void setEdgeVertexFunctor(FunctorType* f) { edgeVertexFunctor = f ; }
	void setFaceVertexFunctor(FunctorType* f) { faceVertexFunctor = f ; }
	void setVolumeVertexFunctor(FunctorType* f) { volumeVertexFunctor = f ; }

	void computeVertexVertexFunctor(Dart d) { (*vertexVertexFunctor)(d); }
	void computeEdgeVertexFunctor(Dart d) { (*edgeVertexFunctor)(d); }
	void computeFaceVertexFunctor(Dart d) { (*faceVertexFunctor)(d); }
	void computerVolumeVertexFunctor(Dart d) { (*volumeVertexFunctor)(d); }

	/*! @name Levels Management
	 *  Operations to manage the levels of an Implicit Hierarchical 3-map
	 *************************************************************************/

	void incCurrentLevel();

	void decCurrentLevel();


	//@{
	//!
	/*!
	 *
	 */
	unsigned int getCurrentLevel() ;

	//!
	/*!
	 *
	 */
	void setCurrentLevel(unsigned int l) ;

	//!
	/*!
	 *
	 */
	unsigned int getMaxLevel() ;

	//!
	/*!
	 *
	 */
	unsigned int getDartLevel(Dart d) ;

	//!
	/*!
	 *
	 */
	void setDartLevel(Dart d, unsigned int i) ;
	//@}

	/*! @name Id Management
	 * Operations to manage the ids of edges and faces
	 *************************************************************************/

	//@{
	//! Give a new unique id to all the edges of the map
	/*!
	 */
	void initEdgeId() ;

	//! Return the next available edge id
	/*!
	 */
	unsigned int getNewEdgeId() ;

	//! Return the id of the edge of d
	/*!
	 */
	unsigned int getEdgeId(Dart d) ;

	//! Set an edge id to all darts from an orbit of d
	/*!
	 */
	void setEdgeId(Dart d, unsigned int i, unsigned int orbit); //TODO a virer
	void setEdgeId(Dart d, unsigned int i);

	//! Give a new unique id to all the faces of the map
	/*!
	 */
	void initFaceId() ;

	//! Return the next available face id
	/*!
	 */
	unsigned int getNewFaceId() ;

	//! Return the id of the face of d
	/*!
	 */
	unsigned int getFaceId(Dart d) ;

	//! Set a face id to all darts from an orbit of d
	/*!
	 */
	void setFaceId(Dart d, unsigned int i, unsigned int orbit); //TODO a virer
	void setFaceId(unsigned int orbit, Dart d);
	//@}

	/*! @name Cells Information
	 * Operations to manage the cells informations :
	 *************************************************************************/

	//@{
	//! Return the level of insertion of the vertex of d
	/*!
	 */
	unsigned int vertexInsertionLevel(Dart d) ;

	//! Return the level of the edge of d in the current level map
	/*!
	 */
	unsigned int edgeLevel(Dart d) ;

	//! Return the level of the face of d in the current level map
	/*!
	 */
	unsigned int faceLevel(Dart d);

	//! Return the level of the volume of d in the current level map
	/*!
	 */
	unsigned int volumeLevel(Dart d);

	//! Return the oldest dart of the face of d in the current level map
	/*!
	 */
	Dart faceOldestDart(Dart d);

	//! Return the oldest dart of the volume of d in the current level map
	/*!
	 */
	Dart volumeOldestDart(Dart d);

	//! Return true if the edge of d in the current level map
	//! has already been subdivided to the next level
	/*!
	 */
	bool edgeIsSubdivided(Dart d) ;

	//! Return true if the edge of d in the current level map
	//! is subdivided to the next level,
	//! none of its resulting edges is in turn subdivided to the next level
	//! and the middle vertex is of degree 2
	/*!
	 */
	bool edgeCanBeCoarsened(Dart d);

	//! Return true if the face of d in the current level map
	//! has already been subdivided to the next level
	/*!
	 */
	bool faceIsSubdivided(Dart d) ;

	//!
	/*!
	 */
	bool faceCanBeCoarsened(Dart d);

	//! Return true if the volume of d in the current level map
	//! has already been subdivided to the next level
	/*!
	 */
	bool volumeIsSubdivided(Dart d);

	//!
	/*!
	 */
	bool volumeIsSubdividedOnce(Dart d);


	/*! @name
	 *************************************************************************/

	//!
	/*!
	 */
	bool neighborhoodLevelDiffersMoreThanOne(Dart d);

	//! wired !!!
	/*!
	 */
	bool coarsenNeighborhoodLevelDiffersMoreThanOne(Dart d);
	//@}

	/*! @name Cell Functors
	 * Redefition of the 3-maps map traversor
	 *************************************************************************/

	//@{
	virtual Dart begin() const;

	virtual Dart end() const;

	virtual void next(Dart& d) const ;

	virtual bool foreach_dart_of_vertex(Dart d, FunctorType& f, unsigned int thread = 0) ;

	virtual bool foreach_dart_of_edge(Dart d, FunctorType& f, unsigned int thread = 0) ;

	bool foreach_dart_of_oriented_face(Dart d, FunctorType& f, unsigned int thread = 0);
	virtual bool foreach_dart_of_face(Dart d, FunctorType& f, unsigned int thread = 0) ;

	bool foreach_dart_of_oriented_volume(Dart d, FunctorType& f, unsigned int thread = 0);
	virtual bool foreach_dart_of_volume(Dart d, FunctorType& f, unsigned int thread = 0) ;

	virtual bool foreach_dart_of_cc(Dart d, FunctorType& f, unsigned int thread = 0) ;


	virtual bool foreach_dart_of_vertex2(Dart d, FunctorType& f, unsigned int thread = 0);

	virtual bool foreach_dart_of_edge2(Dart d, FunctorType& f, unsigned int thread = 0);

	virtual bool foreach_dart_of_face2(Dart d, FunctorType& f, unsigned int thread = 0);
	//@}

	template <unsigned int ORBIT>
		unsigned int getEmbedding(Dart d);
} ;

template <typename T, unsigned int ORBIT>
class AttributeHandler_IHM : public AttributeHandler<T, ORBIT>
{
public:
	typedef T DATA_TYPE ;

	AttributeHandler_IHM() : AttributeHandler<T, ORBIT>()
	{}

	AttributeHandler_IHM(GenericMap* m, AttributeMultiVector<T>* amv) : AttributeHandler<T, ORBIT>(m, amv)
	{}

	AttributeMultiVector<T>* getDataVector() const
	{
		return AttributeHandler<T, ORBIT>::getDataVector() ;
	}

	bool isValid() const
	{
		return AttributeHandler<T, ORBIT>::isValid() ;
	}

	virtual T& operator[](Dart d) ;

	virtual const T& operator[](Dart d) const ;

	T& operator[](unsigned int a)
	{
		return AttributeHandler<T, ORBIT>::operator[](a) ;
	}

	const T& operator[](unsigned int a) const
	{
		return AttributeHandler<T, ORBIT>::operator[](a) ;
	}

} ;

template <typename T>
class VertexAttribute_IHM : public IHM::AttributeHandler_IHM<T, VERTEX>
{
public:
	VertexAttribute_IHM() : IHM::AttributeHandler_IHM<T, VERTEX>() {}
	VertexAttribute_IHM(const IHM::AttributeHandler_IHM<T, VERTEX>& ah) : IHM::AttributeHandler_IHM<T, VERTEX>(ah) {}
	VertexAttribute_IHM<T>& operator=(const IHM::AttributeHandler_IHM<T, VERTEX>& ah) { this->IHM::AttributeHandler_IHM<T, VERTEX>::operator=(ah); return *this; }
};

} //namespace IHM
} // Volume
} //namespace Algo

} //namespace CGoGN

#include "Algo/ImplicitHierarchicalMesh/ihm3.hpp"

#endif

////!
///*!
// */
//bool faceIsSubdividedOnce(Dart d);

