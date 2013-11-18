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

#ifndef __GMAP1_H__
#define __GMAP1_H__

#include "Topology/gmap/gmap0.h"

namespace CGoGN
{

/**
* The class of 1-GMap
*/
class GMap1 : public GMap0
{
protected:
	AttributeMultiVector<Dart>* m_beta1 ;

	void init() ;

public:
	typedef GMap0 ParentMap;

	static const unsigned int DIMENSION = 1 ;



	GMap1();

	virtual std::string mapTypeName() const;

	virtual unsigned int dimension() const;

	virtual void clear(bool removeAttrib);

	virtual void update_topo_shortcuts();

	virtual void compactTopoRelations(const std::vector<unsigned int>& oldnew);

	/*! @name Basic Topological Operators
	 * Access and Modification
	 *************************************************************************/

	virtual Dart newDart();

	Dart beta1(const Dart d);

	template <int N>
	Dart beta(const Dart d);

	Dart phi1(const Dart d);

	Dart phi_1(const Dart d);

	template <int N>
	Dart phi(const Dart d);

	Dart alpha1(const Dart d);

	Dart alpha_1(const Dart d);

protected:
	void beta1sew(Dart d, Dart e);

	void beta1unsew(Dart d);

public:
	/*! @name Constructors and Destructors
	 *  To generate or delete cells in a 1-G-map
	 *************************************************************************/

	//@{
	/**
	* create a new face
	* @param nbEdges the number of sides of face
	* @return a dart of the edge
	*/
	Dart newCycle(unsigned int nbEdges);

	//! Create an new face for boundary (marked)
	/*! @param nbEdges the number of edges
	 *  @return return a dart of the face
	 */
//	Dart newBoundaryCycle(unsigned int nbEdges);

	//! Delete a face erasing all its darts
	/*! @param d a dart of the face
	 */
	void deleteCycle(Dart d) ;
	//@}

	/*! @name Topological Operators
	 *  Topological operations on 1-G-maps
	 *************************************************************************/

	//@{
	//! Cut an edge inserting a new dart between d and its successor in the face
	/*! @param d the edge to cut
	 *  @return a dart of the new vertex
	 * \image hmtl map1_cutEdge.png
	 */
	Dart cutEdge(Dart d);

	//! Undo the cut of the edge of d
	/*! @param d a dart of the edge to uncut
	 */
	void uncutEdge(Dart d);

	//! Collapse an edge of a face
	/*!  \warning Dart d no longer exists after the call
	 *  @param d the edge
	 */
	void collapseEdge(Dart d);

	//! Split a face between vertices d and e
	/*! \pre Dart d and e MUST be different and belong to the same face
	 *  @param d first dart in the face
	 *  @param e second dart in the face
	 */
	void splitCycle(Dart d, Dart e);

	//! Merge the two faces of d and e, darts d & e disappear
	/*! \pre Dart d and e MUST belong to distinct faces
	 *  \warning Darts d and e no longer exist after the call
	 *  @param d a dart in the first face
	 *  @param e a dart in the second face
	 */
	void mergeCycles(Dart d, Dart e);

	//! Link two faces by adding an edge between two vertices
	/*! \pre Dart d and e MUST be different and belong to distinct face
	 *  @param d first dart in the face
	 *  @param e second dart in the face
	 */
	void linkCycles(Dart d, Dart e);
	//@}

	/*! @name Topological Queries
	 *  Return or set various topological information
	 *************************************************************************/

	//@{
	//! Test if dart d and e belong to the same oriented face
	/*! @param d a dart
	 *  @param e a dart
	 */
	bool sameOrientedCycle(Dart d, Dart e) ;

	//! Test if dart d and e belong to the same face
	/*! @param d a dart
	 *  @param e a dart
	 */
	bool sameCycle(Dart d, Dart e) ;

	//! Length of a face (its number of oriented edges)
	/*! @param d a dart of the face
	 *  @return the length of the face
	 */
	unsigned int cycleDegree(Dart d);

	//! Check the Length of a cycle (its number of oriented edges)
	/*! @param d a dart of the cycle
	 *  @param le the length to compare
	 *  @return  negative/null/positive if face degree is less/equal/greater than given degree
	 */
	 int checkCycleDegree(Dart d, unsigned int le) ;

	/**
	 * check if the face of d is a triangle
	 * @return a boolean indicating if the face is a triangle
	 */
	bool isCycleTriangle(Dart d);
	//@}

	/*! @name Cell Functors
	 *  Apply functors to all darts of a cell
	 *************************************************************************/

	//@{
	/**
	* Apply a functor on each dart of a vertex
	* @param d a dart of the vertex
	* @param fonct functor obj ref
	*/
	bool foreach_dart_of_vertex(Dart d, FunctorType& fonct, unsigned int thread=0);

	/**
	* Apply a functor on each dart of an edge
	* @param d a dart of the edge
	* @param fonct functor obj ref
	*/
	bool foreach_dart_of_edge(Dart d, FunctorType& fonct, unsigned int thread=0);

	/**
	* Apply a functor on each dart of an oriented cc (face)
	* @param d a dart of the oriented cc
	* @param fonct functor obj ref
	*/
	bool foreach_dart_of_oriented_cc(Dart d, FunctorType& f, unsigned int thread=0);

	//! Apply a functor on every dart of a cc (face)
	/*! @param d a dart of the cc
	 *  @param f the functor to apply
	 */
	bool foreach_dart_of_cc(Dart d, FunctorType& fonct, unsigned int thread=0);
	//@}
};

} // namespace CGoGN

#include "Topology/gmap/gmap1.hpp"

#endif
