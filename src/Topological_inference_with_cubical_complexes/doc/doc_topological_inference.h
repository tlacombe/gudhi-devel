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
 * There are various objects that introduce intersting topology, and in particular homology and persistent homology. 
 * Typical ones are: sublevelsets of continuous functions, filtered complexes created based 
 * on point clouds, distance functions from sets, distributions and similar. The class Topological_inference_with_cubical_complexes
 * create a cubical complexes with filtration that gives partially constant approximation of such a functions. 
 * In case of very large point clouds in low dimensiona spaces computing persistent homology of a distance from a point cloud 
 * may be superior to use Rips, Cech or Alpha complexes. 
 * 
 * \section Topological_inference_with_cubical_complexes_idea Idea
 * The class Topological_inference_with_cubical_complexes create a filtered cubical complex which provides a partially constant
 * approximation of a function defined on rectangular domain in Euclidean space. To acheive this, first the domain is 
 * covered with a cubical complex. Any implementation of cubical complex that provides certain methods (listed in the refeence 
 * manual) can be used here. Later the approximation of a given function is obtained by computing the value of the function 
 * on the center of every top dimensional cube. 
 * 
 * Here is the list of functions operating on datasets that are implemented as a part of this package:
 * 
 * (1) Function can be given as an algorithm that given a point x, compute f(x). The algorithm can implement a formula, numerical mentod or any other way of assigning f(x) to x.
 * As an example, the figure below illustrates four partially constant approximations of a function assigning to x its distance from a unit circle in the plane. 
 * From left to right, the resolutions of the grid are: 5 by 5, 10 by 10, 50 by 50 and 80 by 80.
 * \image html top_inference.png
 * 
 * (2) Another way is to work with distance functions from a given collection of points. This package implements various distance functions from point clouds. They are 
 * typically build based on a top of a basic (kernel) distance used to compute distance between two points. Currently we support the following kernel distances:  
 * Euclidean, Manhattan and maximum distance.
 * Each of the kernel distances can be computed in a standard Euclidean space or in a periodic domain.
 * The kernel distance is later used to define a function based on a point cloud. Here are the functions of the point clouds that are currently supported
 * (and can be used with any of the kernel distances):
 * 
 * (a) kernels_centerd_in_point_cloud - for any given point x, f(x) is a sum of values of kernel function for x and any other point p from the initial point cloud. 
 * 
 * (b) Sum_of_distances_from_points - for any given point x, f(x) is a sum of distances of x to the points of the point cloud.
 * 
 * (c) Distance_to_k_th_closest_point - for any x, f(x) is a distance of x to the k-th nearest neighbor in the point cloud given a the input. This class uses a brute force quadratic
 * time algorithm, but at the same time allow to use any kernel both on periodic and non periodic domains. 
 * 
 * (d) Distance_to_k_th_closest_point_k_d_tree - for any x, f(x) is a distance of x to the k-th nearest neighbor in the point cloud given a the input. This class require CGAL
 * in version at least 4.8.1 and uses efficient implementation based on k-d-trees. It is hovewer restricted to Euclidean distance in non-periodic domains. 
 * 
 * As an example please consult a periodic function (on a domain [-1.5,1.5]^2) obtained using class Sum_of_distances_from_points
 * where for every point x of a grid f(x) is the sum of Euclidean distances to the point cloud sampled from unit circle.
 * \image html periodic_distance_from_cicrle_nonperiodic_domain.png
 * 
 * (3) Distance function from a collection of cubes that satisfy certain predicate. It can be computed using Morphological_operations_cubical_complex class.
 * This class allows to re-define the filtration of a cubical complex taking into accound the original filtration in the following way. 
 * Initially, once a cubical complex is given or constructed, Morphological_operations_cubical_complex class iterate through all its top dimensional cubes.
 * The ones that satisfy choosen predicate are considered to belong to the 'set'. All other top dimensional cubes became elements of 'set complement'.
 * Given a set and its complement, three basic morphological operations can be performed:
 * (a) Erosion - the elements of the set that are direct neighbors of elements from the set complement are given a new filtration value equal to predefined step_size.
 * The element of the set that did not have already assigned value, and are neighbors of the elements which get the value step_size, will get the value 2*step_size. And so on.
 * (b) Dilation - the elements of the set complement that are direct neighbors of elements from the set are given a new filtration value equal to predefined step_size.
 * The element of the set complement that did not have already assigned value, and are neighbors of the elements which get the value step_size, will get the value 2*step_size. And so on.
 * (c) Both erosion and dilation - In this case, both erosion and dilation is performed at the same time. In this case, the procedure starts from the elements of the set that are 
 * neighbors of the elements of the set complement.
 * 
 * We can choose here between two types of neighborhoods of top dimensional cubes:
 * (a) full_face - neighbors of a top dimensional cube C are all those top dimensional cubes D sharing codimension 1 face with C.
 * (b) all - neighbors of a top dimensional cube C are all those top dimensional cubes D that have nonempty intersection with C.
 * 
 * At the moment the following predicates are implemented:
 * (a) Filtration_above_certain_value - predicate will return true iff given filtration is above certain value. 
 * (b) Filtration_below_certain_value - predicate will return true iff given filtration is below certain value. 
 * (c) Filtration_in_range - predicate will return true iff given filtration is in a given range. 
 * (d) Filtration_equal - predicate will return true iff given filtration is equal certain value. 
 * (e) Always_true  - predicate that will always return true. This is a default tempate parameter of a class Morphological_operations_cubical_complex 
 * Consult Morphological_operations_cubical_complex.h for further details.
 * 
 * An example of the the construction described above is given at the following illustration. On the left, the initial filtration of a two dimensional
 * cubical complex. Second from the left: a result of applying a predicate Filtration_below_certain_value with a parameter 5 to that complex. As a consequence all
 * the cubes with value below 5 became the elements of the set. Third from the left: the result of erosion using full_face neighborhood. Fourth from the left: the result of dilation
 * of the set using full_face neighborhood. In both cases, step_size is set to 1.
 * \image html morphological_op_illustration.png
 * 
 * As an example, let us consider a point cloud obtained in the following way: we sample 50 points from each of the cubes:
 * [-1,0]x[1,2] , [0,1]x[1,2] ,[1,2]x[1,2] 
 * [-1,0]x[0,1] ,             ,[1,2]x[0,1]
 * [-1,0]x[-1,0], [0,1]x[-1,0],[1,2]x[-1,0]
 * and put a 100 by 100 cubical grid to cover the whole [-1,2]x[-1,2] domain. The function giving filtration of this grid is a distance to
 * 5 nearest neighbors in the point cloud. Later we apply Filtration_below_certain_value predicate with parameters (from left to right:
 * 0.1, 0.05, 0.02, 0.015 and 0.01 and dilation starting from those sets. The obtained filtrations are presented in the picture below. 
 * \image html morphological_operation.png
 * 
 * There is another way of constructing objects of Morphological_operations_cubical_complex. One can start with the number of maximal cubes in 
 * the cubical complex and the positions the top dimensional cubes that are supposed to be in the set. The positions are encoded as counters in a way
 * described in the picture below. 
 * \image html counters_of_cubes.png
 * Given such a input, the class Morphological_operations_cubical_complex define elements of sets and its complement. Later the morphological operations
 * can be used for the obtained cubical complexes. 
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
