/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016  INRIA Sophia-Antipolis (France)
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

#ifndef GUDHI_TC_CONFIG_H
#define GUDHI_TC_CONFIG_H

//#include <CGAL/config.h>

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false
// positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

//=========================== Alpha-TC or not? ================================

//#define GUDHI_ALPHA_TC
  //#define GUDHI_USE_A_FIXED_ALPHA
  const double GUDHI_TC_ALPHA_VALUE = 0.3;

//========================= Debugging & profiling =============================
#define GUDHI_TC_PROFILING
#define GUDHI_TC_VERBOSE
//#define GUDHI_TC_VERY_VERBOSE
//#define GUDHI_TC_PERFORM_EXTRA_CHECKS
//#define GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

// Solving inconsistencies: only perturb the vertex, the simplex or more?
#define GUDHI_TC_PERTURB_THE_CENTER_VERTEX_ONLY // Best technique
//#define GUDHI_TC_PERTURB_THE_SIMPLEX_ONLY // Best technique ex-aequo
//#define GUDHI_TC_PERTURB_THE_1_STAR // BAD TECHNIQUE
//#define GUDHI_TC_PERTURB_N_CLOSEST_POINTS // perturb the GUDHI_TC_NUMBER_OF_PERTURBED_POINTS closest points
// Otherwise, perturb one random point of the simplex

// Only used if GUDHI_TC_PERTURB_N_CLOSEST_POINTS is defined
#define GUDHI_TC_NUMBER_OF_PERTURBED_POINTS(intr_dim) (1) // Good technique
//#define GUDHI_TC_NUMBER_OF_PERTURBED_POINTS(intr_dim) (intr_dim + 2) // Quite good technique

#define GUDHI_MESH_D_PROFILING
#define GUDHI_MESH_D_VERBOSE

//========================= Strategy ==========================================
#define GUDHI_TC_NANOFLANN_IS_AVAILABLE
//#define GUDHI_TC_ANN_IS_AVAILABLE
//#define GUDHI_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
//#define GUDHI_TC_PERTURB_EACH_STAR_UNTIL_FIXED
#define GUDHI_TC_GLOBAL_REFRESH
#define GUDHI_TC_PERTURB_POSITION
# define GUDHI_TC_PERTURB_POSITION_TANGENTIAL // default
//# define GUDHI_TC_PERTURB_POSITION_GLOBAL
//#define GUDHI_TC_PERTURB_WEIGHT
//#define GUDHI_TC_PERTURB_TANGENT_SPACE

//========================= Parameters ========================================

// PCA will use BASE_VALUE_FOR_PCA^intrinsic_dim points
const std::size_t BASE_VALUE_FOR_PCA = 5;
const unsigned int BASE_VALUE_FOR_ALPHA_TC_NEIGHBORHOOD = 5;

#endif // GUDHI_TC_CONFIG_H
