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

/**  \defgroup Topological_inference_with_cubical_complexes Topological_inference_with_cubical_complexes distance
 * 
 * \author    Pawel Dlotko
 * @{
 * 
 * \section Topological_inference_with_cubical_complexes_introduction Introduction.
 * 
 * There are various objects that introduce intersting topology and persistent topology. 
 * Typical ones are: sublevelsets of continuous functions, point clouds, distance functions 
 * from sets, distributions and similar. The class Topological_inference_with_cubical_complexes_introduction
 * allows all sorts of those computations. Sometimes using it is superior to use Rips, Cech or Alpha complexes. 
 * 
 * \section Topological_inference_with_cubical_complexes_idea Idea
 * The class Topological_inference_with_cubical_complexes create a cubical complex that provides partially constant
 * approximations of functions defined on rectangular domains in Euclidean spaces. To acheive this, first the domain is 
 * covered by a cubical complex. Any implementation of cubical complex that implements certain methods can be used here.
 * Later the approximation of a given function is obtained by computing the value of the function on the center of every 
 * top dimensional cube. 
 * 
 * There are various ways, some of them listed below, the functions can be obtained. Some of them are implemented in this 
 * package and can be used as a unit. 
 * 
 * (1) Function can be given as an algorithm that given a point x, compute f(x). It can be given by a formula, numerical mentod or any other way.
 * Consult the figure below for a distance function from a unit curcle in the plane. From left to right, the resolutions of the grid are:
 * 5 by 5, 10 by 10, 50 by 50 and 80 by 80.
 * \image html top_inference.png
 * (2) One can consider distance functions from a collection of points. This package implement various distance functions from point cloud. They are 
 * typically build based on a based (kernel) distance. Currently we support the following kernel distances:  Euclidean, Manhattan, maximum distance.
 * Each of the kernel distances can be computed in Euclidean space, or in a periodic domain by using periodic_domain_distance class. 
 * The kernel distance is later used to define a function based on a point cloud. Here are the functions of the point clouds that are currently supported
 * (and can be used with any of the kernel distances):
 * (a) kernels_centerd_in_point_cloud - for any given point p it compute sum of values of kernel function for p and any other point q from the initial point cloud. 
 * (b) Sum_of_distances_from_points - for any given point p it compute sum of distances of p to the points of the point cloud.
 * (c) Distance_to_k_th_closest_point - return the distance to the k-th nearest neighbor in the point cloud given a the input. 
 * As an example please consult the 
 *
 * \image html perturb_pd.png On this picture, the red edges represent the matching. The bottleneck distance is the length of the longest edge.
 *
 */
/** @} */  // end defgroup bottleneck_distance

}  // Topological_inference_with_cubical_complexes persistence_diagram

}  // namespace Gudhi

#endif  // DOC_TOPOLOGICAL_INFERENCE_H_
