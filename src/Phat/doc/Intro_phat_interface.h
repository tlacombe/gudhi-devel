/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DOC_PHAT_INTRO_PHAT_INTERFACE_H_
#define DOC_PHAT_INTRO_PHAT_INTERFACE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace phat_interface {

/** \defgroup phat_interface PHAT interface

  \author    Pawel Dlotko

 \section phatconcept Concept

 We define the concept FilteredComplex which enumerates the requirements for a class to represent a filtered complex
 from which PHAT persistence may be computed.
 We use the vocabulary of simplicial complexes, but the concept is valid for any type of cell complex.
 The main requirements are the definition of:
 \li type <CODE>Simplex_handle</CODE> to manipulate simplices,
 \li type <CODE>Simplex_key</CODE> to manipulate simplices,
 \li type <CODE>Filtration_value</CODE> to manipulate filtration values of simplices,
 \li method <CODE>int num_simplices()</CODE> returning the number of simplices in the simplex,
 \li method <CODE>int dimension()</CODE> returning the simplex dimension,
 \li method <CODE>int dimension(Simplex_handle)</CODE> returning the dimension of a simplex,
 \li method <CODE>assign_key(Simplex_handle, Simplex_key)</CODE> that returns the key associated to a Simplex_handle,
 \li method <CODE>Simplex_key key(Simplex_handle)</CODE> to assign a key to a given Simplex_handle,
 \li method <CODE>Simplex_handle simplex (Simplex_key)</CODE> that returns the Simplex_handle from its Simplex_key,
 \li type and method <CODE>Boundary_simplex_range boundary_simplex_range(Simplex_handle)</CODE> that returns a range
 giving access to the codimension 1 subsimplices of the input simplex, as-well-as the coefficients \f$(-1)^i\f$ in the
 definition of the operator \f$\partial\f$. The iterators have value type <CODE>Simplex_handle</CODE>,
 \li type and method <CODE>Filtration_simplex_range filtration_simplex_range ()</CODE> that returns a range giving
 access to all the simplices of the complex read in the order assigned by the indexing scheme,
 \li type and method <CODE>Filtration_value filtration (Simplex_handle)</CODE> that returns the value of
 the filtration on the simplex represented by the handle.

\section Examples
\subsection cubicalphatpersistence Cubical complex persistence with PHAT interface

\include Phat/cubical_complex_persistence_with_phat.cpp

\subsection ripsphatpersistence Rips complex persistence with PHAT interface

\include Phat/rips_persistence_with_phat.cpp

 \copyright GNU General Public License v3.
 */
/** @} */  // end defgroup phat_interface

}  // namespace phat_interface

}  // namespace Gudhi

#endif  // DOC_PHAT_INTRO_PHAT_INTERFACE_H_
