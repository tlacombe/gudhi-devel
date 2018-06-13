/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

#include <iostream>
#include <fstream>
#include <chrono>
#include <gudhi/Simplex_tree.h>
#include "gudhi/reader_utils.h"

using namespace Gudhi;

typedef Simplex_tree<Simplex_tree_options_zigzag_persistence>  Complex_ds;
typedef Zigzag_edge<Complex_ds>                                Zz_edge_t;

int main()
{
  int dim_max = 3;
  std::vector< Zz_edge_t > edge_filtration;
  //insert manually some edges
  edge_filtration.emplace_back(1,1,   0,true);  //insert {1}
  edge_filtration.emplace_back(2,2,   0,true);  //insert {2}
  edge_filtration.emplace_back(0,0,   0,true);  //insert {0}
  edge_filtration.emplace_back(0,1, 0.1,true);  //insert {0,1}
  edge_filtration.emplace_back(1,2, 0.2,true);  //insert {1,2}
  edge_filtration.emplace_back(3,3,   0,true);  //insert {3}
  edge_filtration.emplace_back(3,0, 0.3,true);  //insert {0,3}
  edge_filtration.emplace_back(0,2, 0.4,true);  //insert {0,2}
  edge_filtration.emplace_back(2,3, 0.5,true);  //insert {2,3}
  edge_filtration.emplace_back(1,3, 0.6,true);  //insert {1,3}
  edge_filtration.emplace_back(0,2,  -1,false); //remove {0,2}
  edge_filtration.emplace_back(1,3,  -1,false); //remove {1,3}
  edge_filtration.emplace_back(0,2,  -1,true);  //insert {0,2}
  edge_filtration.emplace_back(2,3,  -1,false); //remove {2,3}


  std::cout << "Edge filtration: \n";
  for(auto edg : edge_filtration) 
  {
    if(edg.type()) { std::cout << "+ "; } else { std::cout << "- "; }
    std::cout <<  " " << edg.u() << " " << edg.v() << "          " << edg.fil()        << std::endl;
  }
  std::cout << std::endl;

  Complex_ds st;

  // traverse the filtration
  auto zz_rg = st.zigzag_simplex_range(edge_filtration, dim_max);
  
  std::cout << "Simplex filtration: \n";

  for(auto it = zz_rg.begin(); it != zz_rg.end(); ++it ) {
    if(it.arrow_direction()) {std::cout << "+ ";} else {std::cout << "- ";}
    for(auto u : st.simplex_vertex_range(*it)) { std::cout << u << " "; }
      std::cout << "    " << st.filtration(*it) << "\n";
  }
  std::cout << std::endl;
  return 0;
}

