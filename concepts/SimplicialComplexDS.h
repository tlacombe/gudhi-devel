/*
 *  SimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Cl??ment Maria on 12/10/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

/**
 * \brief Data structure for representing a simplicial complex.
 */
struct SimplicialComplexDS
{
	
	/*************************************************/		
	/// \name Objects
	/// @{
	/**
	 * Vertices of a simplicial complex \f$ K = (V,S) \f$
	 */
	typedef unspecified Vertex;
	/**
	 * Simplex type
	 */
	typedef unspecified Simplex;
	///@}
	/*************************************************/	
	
	/*************************************************/		
	/// \name Simplex extra data
	/// 
	/// Extra space allocated for each simplex to store
	/// additional information. If not `void`, it must be persistent,
	/// i.e. it cannot be part of a temporary Simplex object.
	///
	/// @{
	/**
	 * Type of the extra information provided for each simplex
	 * to store additional information
	 */
	typedef unspecified Simplex_data;
	/**
	 * Returns a pointer to the extra information provided for
	 * the simplex s
	 */
	Simplex_data *simplex_data_pointer(Simplex s);
	/// @}
	/*************************************************/	

	
	
	
	/*************************************************/	
	/// \name Iterators
	/// @{
	/**
	 *	Iterator over all simplices of a complex
	 *
	 *	 `value_type` must be a `Simplex`
	 */
	typedef unspecified Complex_simplex_iterator;
	/**
	 *	Iterator over all vertices of a simplex
	 *
	 *	 `value_type` must be a `Vertex`
	 */
	typedef unspecified Simplex_vertex_iterator;
	/**
	 *	 Iterator over all simplices of the boundary of a simplex
	 *
	 *	 `value_type` must be a `Simplex`
	 */
	typedef unspecified Boundary_simplex_iterator;
	/** OPTIONAL
	 * Iterator over all simplices of the coboundary of a simplex
	 *
	 *	`value_type` must be a `Simplex`
	 */
	typedef unspecified Coboundary_simplex_iterator;
	
	
	
	
	/**
	 Returns an iterator to the beginning of the sequence of
	 simplices of a complex
	 */
	Complex_simplex_iterator complex_simplex_iterator_begin(); 
	/**
	 Returns an iterator to the end of the sequence of 
	 simplices of a complex
	 */
	Complex_simplex_iterator complex_simplex_iterator_end(); 
	
	/**
	 * Returns an iterator to the beginning of the sequence of 
	 * vertices of a simplex
	 *
	 * @param s Simplex from which we enumerate the vertices
	 * @return Iterator to the beginning of the sequence of vertices of
	 *	s
	 */
	Simplex_vertex_iterator simplex_vertex_iterator_begin(Simplex s);
	/**
	 Returns an iterator to the end of the sequence of 
	 vertices of a simplex
	 */
	Simplex_vertex_iterator simplex_vertex_iterator_end(Simplex s);
	
	/**
	 Returns an iterator to the beginning of the sequence of
	 simplices of the boundary of a simplex
	 */
	Boundary_simplex_iterator boundary_simplex_iterator_begin(Simplex s);
	/**
	 Returns an iterator to the end of the sequence of 
	 simplices of the boundary of a simplex
	 */
	Boundary_simplex_iterator boundary_simplex_iterator_end(Simplex s);
	
	
	/**
	 Returns an iterator to the beginning of the sequence of
	 simplices of the coboundary of a simplex
	 */
	Coboundary_simplex_iterator coboundary_simplex_iterator_begin(Simplex s);
	/**
	 Returns an iterator to the end of the sequence of 
	 simplices of the coboundary of a simplex
	 */
	Coboundary_simplex_iterator coboundary_simplex_iterator_end(Simplex s);
	
	/// @}
	/*************************************************/	
	
	/*************************************************/	
	/// \name Queries
	/// @{
	/**
	 * Returns the dimension of a simplex s
	 */
	int simplex_dimension(Simplex s);
	/**
	 * Returns the dimension of the complex
	 */
	int complex_dimension();
	
	/// @}
	/*************************************************/	

	/*************************************************/		
	/// \name Modifiers
	/// @{
	/** OPTIONAL
	 * Build the flag complex of dimension maximal_dimension
	 * induced by the 1-skeleton of the simplicial complex
	 */
	void expand(int maximal_dimension);
	/** OPTIONAL
	 * Remove the simplex s and all its cofaces
	 */
	void remove_simplex(Simplex s);
	/** OPTIONAL
	 * Contract an edge s
	 */
	void edge_contraction(Simplex s);
	/// @}
	/*************************************************/	

	

};
