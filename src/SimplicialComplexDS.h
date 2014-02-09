/*
 *  SimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Clément Maria on 12/10/13.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

/**
 * \brief Data structure for representing a simplicial complex.
 *
 * This concept is general, and defines the types Simplex
 * and Vertex.
 */
struct SimplicialComplexDS
{
	/** \brief Defines a type, model of concept MetricSpace.*/
  typedef unspecified MetricSpace;
	/**
	 * \brief Defines a type Vertex, which matches with the one
   * defined in MetricSpace.

	 * The Vertices must admit a total order.
	 */
	typedef MetricSpace::Vertex 		Vertex;
	
	/**
	 * \brief Simplex handle.
	 *
	 * A Simplex_handle represents a unique simplex in the
	 * simplicial complex.
	 */
	typedef unspecified Simplex_handle;
	
	/**
	 * \brief Iterator on the sequence of vertices of
	 * a simplex.
	 *
	 * 'value type' Vertex
	 */
	typedef unspecified Simplex_vertex_iterator;
	
	/** \brief Simplex range type.
	 * 
	 * A range over the vertices of a simplex. A Simplex does not 
	 * need to be a simplex
	 * of the simplicial complex represented. The simplex admits
	 * a canonical orientation induced by the order on its Vertices
	 * in the range returned by simplex_vertex_range(Simplex_handle s)
	 */
	typedef unspecified Simplex_vertex_range;
	
	/**
	 * \brief Returns a range over the Vertices of the input Simplex.
	 *
	 * Specifies the orientation of the simplex corresponding to
	 * the input Simplex_handle.
	 */
	Simplex_vertex_range simplex_vertex_range(Simplex_handle sh);
	///@}
	/*************************************************/	

	/*************************************************/	
	/// \name Queries
	/// @{
	/**
	 * Answer true if the Simplex belongs to the simplicial complex;
	 * false otherwise.
	 *
	 * The set of simplices for which this query answers true defines
	 * the set of simplices of the simplicial complex. It must satisfy 
	 * the subsimplex closeness of simplicial complexes.
	 */
	bool does_simplex_belong_to_complex(Simplex_handle s);
	/**
	 * Returns the dimension of a simplex s, i.e. the
	 * number of Vertices minus 1
	 */
	int simplex_dimension(Simplex_handle s);
	/**
	 * Returns the dimension of the complex, i.e. the maximal
	 * dimension of a simplex in the simplicial complex.
	 */
	int complex_dimension();
	/// @}
	/*************************************************/	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 *	Iterator over all simplices of a complex
	 *
	 *	 `value_type` must be a `Simplex`
	 */
//	typedef unspecified Complex_simplex_iterator;
	/**
	 *	Range of all simplices of a complex, corresponding to the
	 *	iterator type `Complex_simplex_iterator`
	 */
//	typedef unspecified Complex_simplex_range;
	/** OPTIONAL
	 * Iterator over all simplices of the coboundary of a simplex
	 *
	 *	`value_type` must be a `Simplex`
	 */
//	typedef unspecified Coboundary_simplex_iterator;
	/// @}
	
	/*************************************************/	

	/**
	 Returns an iterator to the beginning of the sequence of
	 simplices of a complex
	 */
//	Complex_simplex_iterator complex_simplex_begin(); 
	/**
	 Returns an iterator to the end of the sequence of 
	 simplices of a complex
	 */
//	Complex_simplex_iterator complex_simplex_end(); 
	/**
	 Returns a range of the sequence of simplices of a complex
	 */
//	Complex_simplex_range complex_simplices(); 

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
	

	/*************************************************/		
	/// \name Modifiers
	/// @{
	/** OPTIONAL
	 * Build the flag complex of dimension maximal_dimension
	 * induced by the 1-skeleton of the simplicial complex
	 */
//	void expand(int maximal_dimension);
	/** OPTIONAL
	 * Remove the simplex s and all its cofaces
	 */
//	void remove_simplex(Simplex s);
	/** OPTIONAL
	 * Contract an edge s
	 */
//	void edge_contraction(Simplex s);
	/// @}
	/*************************************************/	

	/// \name Simplex extra data accessor
	/// @{
	/**
	 * Returns a pointer to the extra information provided for
	 * the simplex s.  If not `Simplex_data` is not `void`, it must be
	 * persistent, i.e. it cannot be part of a temporary Simplex object.
	 */
//	Simplex_data *simplex_data_pointer(Simplex s);
	/// @}
};
