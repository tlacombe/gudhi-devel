/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2017  INRIA (France)
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

#define BOOST_PARAMETER_MAX_ARITY 12

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Persistent_cohomology.h>

#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>

using Distance_table = std::vector<std::vector<double>>;
Distance_table dist_table;

void initialize_table()
{
  // Example with 5 points with a given distance
  std::vector<double> line0 = {0, 1, 2, 3, 4}; dist_table.push_back(line0);
  std::vector<double> line1 = {1, 0, 1, 2, 3}; dist_table.push_back(line1);
  std::vector<double> line2 = {2, 1, 0, 1, 2}; dist_table.push_back(line2);
  std::vector<double> line3 = {3, 2, 1, 0, 1}; dist_table.push_back(line3);
  std::vector<double> line4 = {4, 3, 2, 1, 0}; dist_table.push_back(line4);
}

double distance(int p1, int p2)
{
  return dist_table[p1][p2];
}


template <typename Distance_table,
          typename Nearest_landmark_table>
void construct_nearest_landmark_table(Distance_table const &dt_,
                                      Nearest_landmark_table &nearest_neighbors) {
  
  typedef std::pair<int, double> id_dist_pair;
  typedef bool (*comp)(id_dist_pair, id_dist_pair);
  for (auto dt_it = dt_.begin(); dt_it != dt_.end(); ++dt_it) {
    // for (int points_i = 0; points_i < nbP; points_i++) {
    typename Nearest_landmark_table::value_type nearest_landmark_range;
    std::priority_queue<id_dist_pair, std::vector<id_dist_pair>, comp> l_heap([](id_dist_pair j1, id_dist_pair j2) {
        return j1.second > j2.second;
      });
    auto landmarks_it = dt_it->begin();
    int landmarks_i = 0;
    for (landmarks_i = 0; landmarks_it != dt_it->end();
         ++landmarks_it, landmarks_i++) {
      id_dist_pair id_dist = std::make_pair(landmarks_i, *landmarks_it);
      l_heap.push(id_dist);
    }
    while (!l_heap.empty()) {
      id_dist_pair id_dist = l_heap.top();
      nearest_landmark_range.push_back(id_dist);
      l_heap.pop();
    }
    nearest_neighbors.push_back(nearest_landmark_range);
  }
}

int main(int argc, char * const argv[]) {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;
  using Witness_complex = Gudhi::witness_complex::Witness_complex<Nearest_landmark_table>;
  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

  Simplex_tree simplex_tree;
  Nearest_landmark_table nlt;

  initialize_table();
  construct_nearest_landmark_table(dist_table, nlt);
  
  Witness_complex witness_complex(nlt);
  witness_complex.create_complex(simplex_tree, 4.1);

  std::cout << "Number of simplices: " << simplex_tree.num_simplices() << std::endl;

  Persistent_cohomology pcoh(simplex_tree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(11);

  pcoh.compute_persistent_cohomology(-0.1);
  pcoh.output_diagram();
}
