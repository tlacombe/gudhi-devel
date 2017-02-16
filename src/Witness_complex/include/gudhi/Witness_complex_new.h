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

#ifndef WITNESS_COMPLEX_NEW_H_
#define WITNESS_COMPLEX_NEW_H_

#include <gudhi/Witness_complex.h>
#include <gudhi/Kd_tree_search.h>
#include <gudhi/Active_witness/Active_witness.h>
#include <gudhi/Active_witness/Witness_for_simplex.h>
#include <gudhi/Active_witness/Sib_vertex_pair.h>
#include <gudhi/Witness_complex/all_faces_in.h>
#include <gudhi/Witness_complex/check_if_neighbors.h>
#include <gudhi/Simplex_tree/Vertex_subtree_iterator.h>
#include <gudhi/Simplex_tree/Fixed_dimension_iterator.h>

#include <utility>
#include <vector>
#include <list>
#include <map>
#include <limits>

namespace Gudhi {
  
namespace witness_complex {

/**
 * \private
 * \class Witness_complex
 * \brief Constructs (weak) witness complex for a given table of nearest landmarks with respect to witnesses.
 * \ingroup witness_complex
 *
 * \tparam Nearest_landmark_table_ needs to be a range of a range of pairs of nearest landmarks and distances.
 *         The range of pairs must admit a member type 'iterator'. The dereference type 
 *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'.
*/
template< class Nearest_landmark_table_ >
class Witness_complex_new : public Witness_complex<Nearest_landmark_table_> {
private:
  typedef typename Nearest_landmark_table_::value_type               Nearest_landmark_range;
  typedef std::size_t                                                Witness_id;
  typedef std::size_t                                                Landmark_id;
  typedef std::pair<Landmark_id, double>                             Id_distance_pair;
  typedef Active_witness<Id_distance_pair, Nearest_landmark_range>   ActiveWitness;
  typedef std::list< ActiveWitness >                                 ActiveWitnessList;
  typedef std::vector< Landmark_id >                                 typeVectorVertex;
  typedef std::vector<Nearest_landmark_range>                        Nearest_landmark_table_internal;

  typedef Landmark_id Vertex_handle;

 protected:
  Nearest_landmark_table_internal              nearest_landmark_table_;
  
 public:
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */

  //@{

  Witness_complex_new()
  {
  }
  
  /**
   *  \brief Initializes member variables before constructing simplicial complex.
   *  \details Records nearest landmark table.
   *  @param[in] nearest_landmark_table needs to be a range of a range of pairs of nearest landmarks and distances.
 *         The range of pairs must admit a member type 'iterator'. The dereference type 
 *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'.
   */

  Witness_complex_new(Nearest_landmark_table_ const & nearest_landmark_table)
    : nearest_landmark_table_(std::begin(nearest_landmark_table), std::end(nearest_landmark_table))
  {    
  }

    
  /** \brief Outputs the (weak) witness complex of relaxation 'max_alpha_square'
   *         in a simplicial complex data structure.
   *  \details The function returns true if the construction is successful and false otherwise.
   *  @param[out] complex Simplicial complex data structure compatible which is a model of
   *              SimplicialComplexForWitness concept.
   *  @param[in] max_alpha_square Maximal squared relaxation parameter.
   *  @param[in] limit_dimension Represents the maximal dimension of the simplicial complex
   *         (default value = no limit).
   */
  template < typename SimplicialComplexForWitness >
  bool create_complex(SimplicialComplexForWitness& complex,
                      double  max_alpha_square,
                      std::size_t limit_dimension = std::numeric_limits<std::size_t>::max()) const      
  {
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
    typedef std::vector<Vertex_handle> Vertex_vector;

    typedef Gudhi::Simplex_tree_fixed_dimension_iterator<SimplicialComplexForWitness> Fixed_dimension_iterator;
    typedef boost::iterator_range<Fixed_dimension_iterator> Fixed_dimension_range;

    typedef Witness_for_simplex<typename ActiveWitness::iterator, ActiveWitnessList> Witnessed_simplex;
    typedef std::list<Witnessed_simplex> Witnessed_simplex_list;
    typedef std::map<Simplex_key, Witnessed_simplex_list> Simplex_witness_list_map;
    
    if (complex.num_vertices() > 0) {
      std::cerr << "Witness complex cannot create complex - complex is not empty.\n";
      return false;
    }
    if (max_alpha_square < 0) {
      std::cerr << "Witness complex cannot create complex - squared relaxation parameter must be non-negative.\n";
      return false;
    }
    if (limit_dimension < 0) {
      std::cerr << "Witness complex cannot create complex - limit dimension must be non-negative.\n";
      return false;
    }
    ActiveWitnessList active_witnesses;
    for (auto w: nearest_landmark_table_)
      active_witnesses.push_back(ActiveWitness(w));

    Simplex_witness_list_map* prev_dim_map = new Simplex_witness_list_map();
    fill_vertices(max_alpha_square, complex, active_witnesses, prev_dim_map);
    
    std::cout << "Size of the vertex map: " << prev_dim_map->size() << std::endl;
    //----------------------- TEST
    for (auto vw_pair: *prev_dim_map) {
      std::cout << "*";
      for (auto v: complex.simplex_vertex_range(vw_pair.first.simplex_handle()))
        std::cout << v;
      std::cout << "\n[";
      for (auto w: vw_pair.second) {
        std::cout << "(" << &(*w.witness_)
                  << ", " << w.last_it_->first
                  << ", " << w.limit_distance_
                  << ")";
      }
      std::cout << "]\n";
    }
    //---------------------------

    Simplex_witness_list_map* dim1_map = new Simplex_witness_list_map();
    fill_edges(max_alpha_square, complex, active_witnesses, prev_dim_map, dim1_map);
    delete prev_dim_map;
    prev_dim_map = dim1_map;

    //----------------------- TEST
    for (auto vw_pair: *prev_dim_map) {
      std::cout << "*";
      for (auto v: complex.simplex_vertex_range(vw_pair.first.simplex_handle()))
        std::cout << v;
      std::cout << "\n[";
      for (auto w: vw_pair.second) {
        std::cout << "(" << &(*w.witness_)
                  << ", " << w.last_it_->first
                  << ", " << w.limit_distance_
                  << ")";        
      }
      std::cout << "]\n";
    }

    Landmark_id k = 2; /* current dimension in iterative construction */
    while (!active_witnesses.empty() && k <= limit_dimension) {
      Simplex_witness_list_map* curr_dim_map = new Simplex_witness_list_map();
      fill_simplices(max_alpha_square, k, complex, active_witnesses, prev_dim_map, curr_dim_map);
      delete prev_dim_map;
      prev_dim_map = curr_dim_map;
      //----------------------- TEST
      for (auto vw_pair: *prev_dim_map) {
        std::cout << "*";
        for (auto v: complex.simplex_vertex_range(vw_pair.first.simplex_handle()))
          std::cout << v;
        std::cout << "\n[";
        for (auto w: vw_pair.second) {
          std::cout << "(" << &(*w.witness_)
                    << ", " << w.last_it_->first
                    << ", " << w.limit_distance_
                    << ")";        
        }
        std::cout << "]\n";
        //--------------------------------
      }
      std::cout << "Active witness list size = " << active_witnesses.size() << "\n";
      k++;
    }
    delete prev_dim_map;
    complex.set_dimension(k-1);
    return true;
  }

  //@}

 private:

  /* \brief Fills the map "vertex -> witnesses for vertices"
   * It is necessary to go through the aw list two times,
   * because of the use of Simplex_handle as key in the map.
   * With every insertion of a vertex, the previous Simplex_handles
   * are invalidated, therefore it is not possible to build 
   * Simplex_tree and the map at the same time.
   */
  template < typename SimplicialComplexForWitness,
             typename ActiveWitnessList,
             typename SimplexWitnessMap >
  void fill_vertices(const double alpha2,
                     SimplicialComplexForWitness& sc,
                     ActiveWitnessList& aw_list,
                     SimplexWitnessMap* sw_map) const
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
    typedef std::vector<Vertex_handle> Vertex_vector;
    typedef typename SimplexWitnessMap::mapped_type Simplex_witness_list;
    typedef typename Simplex_witness_list::value_type WitnessForSimplex;

    for (auto aw_it = aw_list.begin(); aw_it != aw_list.end(); ++aw_it) {
      typename ActiveWitness::iterator l_it = aw_it->begin();
      typename ActiveWitness::iterator end = aw_it->end();
      double filtration_value = 0;
      double norelax_dist2 = std::numeric_limits<double>::infinity();
      for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
        if (l_it->second > norelax_dist2)
          filtration_value = l_it->second - norelax_dist2;
        else
          norelax_dist2 = l_it->second;
        sc.insert_simplex(Vertex_vector(1, l_it->first), filtration_value);
      }
    }
    for (auto aw_it = aw_list.begin(); aw_it != aw_list.end(); ++aw_it) {
      typename ActiveWitness::iterator l_it = aw_it->begin();
      typename ActiveWitness::iterator end = aw_it->end();
      double norelax_dist2 = std::numeric_limits<double>::infinity();
      for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
        Simplex_handle sh = sc.find(Vertex_vector(1, l_it->first));
        Siblings* sib = sc.self_siblings(sh);
        Vertex_handle v = sh->first;
        Simplex_key sk(sib,v);
        (*sw_map)[sk].emplace_back(WitnessForSimplex(l_it, aw_it, norelax_dist2));
        aw_it->increase();
        if (l_it->second < norelax_dist2)
          norelax_dist2 = l_it->second;
      }
    }
  }

  /* \brief Fills the map "edges -> witnesses for edges"
   */
  template < typename SimplicialComplexForWitness,
             typename ActiveWitnessList,
             typename SimplexWitnessMap >
  void fill_edges(const double alpha2,
                  SimplicialComplexForWitness& sc,
                  ActiveWitnessList& aw_list,
                  SimplexWitnessMap* dim0_map,
                  SimplexWitnessMap* dim1_map) const
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
    typedef std::vector<Vertex_handle> Vertex_vector;
    typedef typename SimplexWitnessMap::mapped_type Simplex_witness_list;
    typedef typename Simplex_witness_list::value_type WitnessForSimplex;

    for (auto vw_pair: *dim0_map) {
      for (auto w: vw_pair.second) {
        typename ActiveWitness::iterator l_it = w.last_it_;
        typename ActiveWitness::iterator end = w.witness_->end();
        if (l_it != end)
          l_it++;
        double filtration_value = 0;
        double norelax_dist2 = w.limit_distance_;
        for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
          if (l_it->second > norelax_dist2)
            filtration_value = l_it->second - norelax_dist2;
          else
            norelax_dist2 = l_it->second;
          auto sv_range = sc.simplex_vertex_range(vw_pair.first.simplex_handle());
          Vertex_vector vertices(sv_range.begin(), sv_range.end());
          vertices.push_back(l_it->first);
          sc.insert_simplex(vertices, filtration_value);
        }
      }
    }
    for (auto vw_pair: *dim0_map) {
      for (auto w: vw_pair.second) {
        typename ActiveWitness::iterator l_it = w.last_it_;
        typename ActiveWitness::iterator end = w.witness_->end();
        if (l_it != end)
          l_it++;
        double norelax_dist2 = w.limit_distance_;
        for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
          auto sv_range = sc.simplex_vertex_range(vw_pair.first.simplex_handle());
          Vertex_vector vertices(sv_range.begin(), sv_range.end());
          vertices.push_back(l_it->first);
          Simplex_handle sh = sc.find(vertices);
          Siblings* sib = sc.self_siblings(sh);
          Vertex_handle v = sh->first;
          Simplex_key sk(sib,v);
          (*dim1_map)[sk].emplace_back(WitnessForSimplex(l_it, w.witness_, norelax_dist2));
          w.witness_->increase();
          if (l_it->second < norelax_dist2)
            norelax_dist2 = l_it->second;
        }
        w.witness_->decrease();
        if (w.witness_->counter() == 0)
          aw_list.erase(w.witness_);
      }
    }
  }

  /* \brief Fills the map "k-simplex -> witnesses for edges"
   */
  template < typename SimplicialComplexForWitness,
             typename ActiveWitnessList,
             typename SimplexWitnessMap >
  void fill_simplices(double alpha2,
                      std::size_t k,
                      SimplicialComplexForWitness& sc,
                      ActiveWitnessList& aw_list,
                      SimplexWitnessMap* prev_dim_map,
                      SimplexWitnessMap* curr_dim_map) const
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
    typedef std::vector<Vertex_handle> Vertex_vector;
    typedef typename SimplexWitnessMap::mapped_type Simplex_witness_list;
    typedef typename Simplex_witness_list::value_type WitnessForSimplex;
    typedef Gudhi::Simplex_tree_vertex_subtree_iterator<SimplicialComplexForWitness> Vertex_subtree_iterator;
    typedef boost::iterator_range<Vertex_subtree_iterator> Vertex_subtree_range;

    
    for (auto sw_pair: *prev_dim_map) {
      //----------------- TEST
      std::cout << "*";
      for (auto v: sc.simplex_vertex_range(sw_pair.first.simplex_handle()))
        std::cout << v;
      std::cout << std::endl;
      //-----------------
      auto v_it = sc.simplex_vertex_range(sw_pair.first.simplex_handle()).begin();
      unsigned counter = 0; /* I need the first vertex before last */ 
      while (counter != k-2) {
        v_it++;
        counter++;
      }
      Vertex_handle v1 = *(v_it++);
      Vertex_handle v0 = *v_it;
      Vertex_subtree_range v0s_range(Vertex_subtree_iterator(&sc, v0, k-1),
                                    Vertex_subtree_iterator());
      Vertex_subtree_range v1s_range(Vertex_subtree_iterator(&sc, v1, k-1),
                                    Vertex_subtree_iterator());
      for (auto sh2: v0s_range) {
        //----------------- TEST
        for (auto v: sc.simplex_vertex_range(sh2))
          std::cout << v;
        std::cout << std::endl;
        //-----------------
        Vertex_vector coface;
        if (check_if_neighbors(sc, sw_pair.first.simplex_handle(), sh2, coface)) {
           //---------------TEST
          std::cout << "Coface: ";
          for (auto v: coface)
            std::cout << v;
          std::cout << std::endl;
          //-----------------
          double filtration_value = 0;
          if (all_faces_in(coface, &filtration_value, sc)) {
            std::pair<Simplex_handle, bool> sh_bool = sc.insert_simplex(coface);
            if (sh_bool.second) {
              Siblings* sib = sc.self_siblings(sh_bool.first);
              Vertex_handle v = sh_bool.first->first;
              curr_dim_map->emplace(Simplex_key(sib,v), Simplex_witness_list());
              std::cout << "Inserted\n";
            }
          }
        }
      }
      for (auto sh2: v1s_range) {
        //----------------- TEST
        for (auto v: sc.simplex_vertex_range(sh2))
          std::cout << v;
        std::cout << std::endl;
        //-----------------
        Vertex_vector coface;
        if (check_if_neighbors(sc, sw_pair.first.simplex_handle(), sh2, coface)) {
          //---------------TEST
          std::cout << "Coface: ";
          for (auto v: coface)
            std::cout << v;
          std::cout << std::endl;
          //-----------------
          double filtration_value = std::numeric_limits<double>::infinity();
          if (all_faces_in(coface, &filtration_value, sc)) {
            std::pair<Simplex_handle, bool> sh_bool = sc.insert_simplex(coface);
            if (sh_bool.second) {
              Siblings* sib = sc.self_siblings(sh_bool.first);
              Vertex_handle v = sh_bool.first->first;
              curr_dim_map->emplace(Simplex_key(sib,v), Simplex_witness_list());
              std::cout << "Inserted\n";
            }
          }
        }
      }
    }
    for (auto vw_pair: *prev_dim_map) {
      for (auto w: vw_pair.second) {
        typename ActiveWitness::iterator l_it = w.last_it_;
        typename ActiveWitness::iterator end = w.witness_->end();
        if (l_it != end)
          l_it++;
        double norelax_dist2 = w.limit_distance_;
        for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
          auto sv_range = sc.simplex_vertex_range(vw_pair.first.simplex_handle());
          Vertex_vector vertices(sv_range.begin(), sv_range.end());
          vertices.push_back(l_it->first);
          Simplex_handle sh = sc.find(vertices);
          if (sh != sc.null_simplex()) {
            Siblings* sib = sc.self_siblings(sh);
            Vertex_handle v = sh->first;
            Simplex_key sk(sib,v);
            (*curr_dim_map)[sk].emplace_back(WitnessForSimplex(l_it, w.witness_, norelax_dist2));
            sc.insert_simplex(vertices, l_it->second - norelax_dist2); // Update the filtration
            w.witness_->increase();
            if (l_it->second < norelax_dist2)
              norelax_dist2 = l_it->second;
          }
        }
        w.witness_->decrease();
        if (w.witness_->counter() == 0)
          aw_list.erase(w.witness_);
      }
    }
    for (auto vm_it = curr_dim_map->begin(); vm_it != curr_dim_map->end(); ++vm_it) {
      if (vm_it->second.empty()) {
        sc.remove_maximal_simplex(vm_it->first.simplex_handle());
        curr_dim_map->erase(vm_it++);
      } 
    }
  }  

};
  
}  // namespace witness_complex

}  // namespace Gudhi

#endif  // WITNESS_COMPLEX_NEW_H_
