/*
 *  SimplexDataFilteredSimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Cl??ment Maria on 12/10/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

/**
 *  \brief Data structure for representing a filtered simplicial complex, where
 *  additional data may be stored in the simplices.
 *
 *  It is designed to be used by any algorithm to compute persistent homology
 *  and persistent cohomology when additional data may be stored in the simplices
 *	by the persistence algorithm.
 *
 *  \extends FilteredSimplicialComplexDS
 */
struct SimplexDataFilteredSimplicialComplexDS
{
	/*************************************************/			
	/// \name Simplex extra data
	/// 
	/// @{
	/**
	 * Type of the extra space provided for each simplex
	 * to store additional information
	 * 
	 * @todo Should we specify methods to access and modify Simplex_data? 
	 */
	typedef unspecified Simplex_data;
	/// @}
	/*************************************************/
};
