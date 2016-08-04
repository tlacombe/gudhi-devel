/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_
#define DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_

// needs namespaces for Doxygen to link on classes
namespace Gudhi {
namespace tangential_complex {

/**  \defgroup tangential_complex Tangential complex
 * 
 * \author    Cl&eacute;ment Jamin
 * 
 * @{
 * 
 * \section definition Definition
 * 
 * A Tangential Delaunay complex is a <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
 * designed to reconstruct a $k$-dimensional manifold embedded in $d$-dimensional Euclidean space. 
 * The input is a point sample coming from an unknown manifold.
 * The running time depends only linearly on the extrinsic dimension $d$
 * and exponentially on the intrinsic dimension $k$.
 *
 * An extensive description of the Tangential complex can be found in \cite tangentialcomplex2014.
 * 
 * \image html "tc_examples.png" "Examples of Tangential complexes"
 * 
 * \section simple_example Simple example
 * 
 * This example builds the Tangential complex of point set.
 * 
 * \include Tangential_complex/example_basic.cpp
 * 
 * \section example_with_perturb Example with perturbation
 * 
 * This example builds the Tangential complex of a point set, then tries to solve inconsistencies
 * by perturbing the positions of points involved in inconsistent simplices.
 * 
 * \include Tangential_complex/example_with_perturb.cpp
 * 
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup tangential_complex

}  // namespace tangential_complex

}  // namespace Gudhi

#endif  // DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_
