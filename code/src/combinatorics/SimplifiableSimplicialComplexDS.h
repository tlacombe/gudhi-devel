/*
 *  SimplifiableSimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Clément Maria on 12/10/13.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

/**
 * \brief Data structure for representing a simplicial complex that
 * can be simplified.
 *
 * The concept provides tests to "control" the changes in topology
 * induced by the simplifications.
 * A model of this concept may be used to compute the homology of the
 * simplicial complex. One would use the Complex_simplex_iterator
 * to convert the simplicial complex into a matrix. I put the 
 * Boundary_simplex_iterator OPTIONAL for now...
 *
 * \extends SimplicialComplexDS
 */
struct SimplifiableSimplicialComplexDS
{
	/*************************************************/		
	/// \name Objects
	/// @{
	/**
	 * Edge type.
	 *
	 */
	typedef unspecified Edge;

	/*************************************************/	
	/// \name Simplicial Complex Iterator
	/// @{
	/**
	 *	Iterator over all simplices of a complex
	 *
	 *	 `value_type` must be a `Simplex`
	 */
	typedef unspecified Complex_simplex_iterator;

	/**
	 * Returns a range of the sequence of
	 * k-simplices of a complex
	 */
	Complex_simplex_iterator simplex_range(int k);


	/**
	 * Returns a range of the sequence of
	 * simplices of a complex
	 */
	Complex_simplex_iterator simplex_range();

	/**
	 * Returns an iterator to the end of the sequence of 
	 * simplices of a complex
	 */
	Complex_simplex_iterator simplex_end();


		/// @}
	/*************************************************/	

	/*************************************************/		
	/// \name OPTIONAL Boundary Iterator
	/// @{
	/**
	 *	 Iterator over all simplices of the boundary of a simplex, i.e.
	 *   the set of codimension $1$ subsimplices of the Simplex.
	 *
	 *	 `value_type` must be a `Simplex`.
	 */
	typedef unspecified Boundary_simplex_iterator;
	/**
	 * Returns an iterator to the beginning of the sequence of
	 * simplices of the boundary of a simplex
	 *
	 * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
	 * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
	 * simplices of the boundary in the order: 
	 * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
	 *
	 * We note that the alternate sum of the simplices given by the iterator
	 * gives the chains corresponding to the boundary of the simplex.
	 */
	Boundary_simplex_iterator boundary_simplex_begin(Simplex s);
	/**
	 * Returns an iterator to the end of the sequence of 
	 * simplices of the boundary of a simplex
	 */
	Boundary_simplex_iterator boundary_simplex_end(Simplex s);
	/// @}
	/*************************************************/	

	/*************************************************/	
	/// \name Simplification
	/// @{
	/**
	 * Remove a Simplex and all its cofaces
	 */
	int remove(Simplex s);

	/**
	 * Check if removing s would modify the homotopy
	 * type of the simplicial complex
	 *
	 * The answer may be yes, no, maybe are other 
	 * (preserve homoemorphism, ...)
	 */
	enum does_removal_preserve_homotopy(Simplex s);


	/**
	 * Proceed to an edge contraction
	 */
	int contract_edge(Edge e);

	/**
	 * Check if contracting e would modify the homotopy
	 * type of the simplicial complex
	 *
	 * The answer may be yes, no, maybe are other 
	 * (preserve homoemorphism, ...)
	 *
	 * @remark We do not talk about the link condition. It's up to
	 * the developper to choose which kind of topological test he wants
	 * to implement, depending on the precision of the test and its
	 * computational cost.
	 */
	enum does_contraction_preserve_homotopy(Edge e);
	/// @}
	/*************************************************/	





	/*************************************************/		
	/// \name OPTIONAL Coboundary Iterator
	/// @{
	/** OPTIONAL
	 * Iterator over all simplices of the coboundary of a simplex
	 *
	 *	`value_type` must be a `Simplex`
	 */
	//	typedef unspecified Coboundary_simplex_iterator;
	/**
	 Returns an iterator to the beginning of the sequence of
	 simplices of the coboundary of a simplex
	 */
	//	Coboundary_simplex_iterator coboundary_simplex_begin(Simplex s);
	/**
	 Returns an iterator to the end of the sequence of 
	 simplices of the coboundary of a simplex
	 */
	//	Coboundary_simplex_iterator coboundary_simplex_end(Simplex s);
	/// @}
	/*************************************************/	
};
