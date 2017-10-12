/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University
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

#ifndef DOC_TOPOLOGICAL_INFERENCE_H_
#define DOC_TOPOLOGICAL_INFERENCE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace Topological_inference_with_cubical_complexes {

/**  \defgroup Topological_inference_with_cubical_complexes Topological_inference_with_cubical_complexes
 * 
 * \author    Pawel Dlotko
 * @{
 * 
 * \section Topological_inference_with_cubical_complexes_introduction Introduction.
 * 
 * There are various objects that introduce intersting homology and persistent homology. 
 * Typical ones are: sublevelsets of continuous functions, filtered complexes created based 
 * in point clouds, distance functions from sets, distributions and similar. The class Topological_inference_with_cubical_complexes_introduction
 * allows create a cubical complexes that can approximate all those functions. 
 * Sometimes using it to approximate distance from a point cloud may be superior to use Rips, Cech or Alpha complexes for very large complexes. 
 * 
 * \section Topological_inference_with_cubical_complexes_idea Idea
 * The class Topological_inference_with_cubical_complexes create a cubical complex that provides partially constant
 * approximations of functions defined on rectangular domains in Euclidean spaces. To acheive this, first the domain is 
 * covered by a cubical complex. Any implementation of cubical complex that implements certain methods can be used here.
 * Later the approximation of a given function is obtained by computing the value of the function on the center of every 
 * top dimensional cube. 
 * 
 * There are various ways, the functions we use to compute values on cubical grid can be obtained. 
 * A lot of them are implemented in this package:
 * 
 * (1) Function can be given as an algorithm that given a point x, compute f(x). The algorithm can implement a formula, numerical mentod or any other way of getting f(x) given x.
 * Consult the figure below for a filtered cubical complexes approximating distance function from a unit circle in the plane. From left to right, the resolutions of the grid are:
 * 5 by 5, 10 by 10, 50 by 50 and 80 by 80.
 * \image html top_inference.png
 * (2) One can consider distance functions from a collection of points. This package implement various distance functions from point cloud. They are 
 * typically build based on a basic (kernel) distance. Currently we support the following kernel distances:  Euclidean, Manhattan and maximum distance.
 * Each of the kernel distances can be computed in a standard Euclidean space or in a periodic domain.
 * The kernel distance is later used to define a function based on a point cloud. Here are the functions of the point clouds that are currently supported
 * (and can be used with any of the kernel distances):
 * 
 * (a) kernels_centerd_in_point_cloud - for any given point p it compute sum of values of kernel function for p and any other point q from the initial point cloud. 
 * 
 * (b) Sum_of_distances_from_points - for any given point p it compute sum of distances of p to the points of the point cloud.
 * 
 * (c) Distance_to_k_th_closest_point - return the distance to the k-th nearest neighbor in the point cloud given a the input. 
 * As an example please consult the periodic (on a domain [-1.5,1.5]^2) function obtained using class Sum_of_distances_from_points
 * where for every point of a grid the sum of Euclidean distances to the point cloud sampled from unit circle have been computed.
 * \image html periodic_distance_from_cicrle_nonperiodic_domain.png
 * 
 * (3) Distance function from a collection of cubes that satisfy certain predicate. TODO
 * 
 * All the cubical complexes used here may, or may not have periodic boundary conditions imposed. Please consult examples and utilities for further details. 
 * 
 *
 *
 */
/** @} */   // end defgroup Topological_inference_with_cubical_complexes

}  // Topological_inference_with_cubical_complexes

}  // namespace Gudhi

#endif  // DOC_TOPOLOGICAL_INFERENCE_H_
