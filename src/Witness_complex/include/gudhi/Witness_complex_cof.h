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

#ifndef WITNESS_COMPLEX_COF_H_
#define WITNESS_COMPLEX_COF_H_

#include <gudhi/Witness_complex.h>
#include <gudhi/Witness_complex_new.h>
#include <gudhi/Kd_tree_search.h>
#include <gudhi/Active_witness/Active_witness.h>
#include <gudhi/Active_witness/Witness_for_simplex.h>
#include <gudhi/Active_witness/Sib_vertex_pair.h>
#include <gudhi/Witness_complex/all_faces_in.h>
#include <gudhi/Witness_complex/check_if_neighbors.h>
#include <gudhi/Simplex_tree/Vertex_subtree_iterator.h>
#include <gudhi/Simplex_tree/Fixed_dimension_iterator.h>

#include <gudhi/Simplex_tree.h>

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
class Witness_complex_cof : public Witness_complex<Nearest_landmark_table_> {
private:
  typedef typename Nearest_landmark_table_::value_type               Nearest_landmark_range;
  typedef std::size_t                                                Witness_id;
  typedef std::size_t                                                Landmark_id;
  typedef std::pair<Landmark_id, double>                             Id_distance_pair;
  typedef Active_witness<Id_distance_pair, Nearest_landmark_range>   ActiveWitness;
  typedef std::list< ActiveWitness >                                 ActiveWitnessList;
  typedef std::vector< Landmark_id >                                 Vertex_vector;
  typedef std::vector<Nearest_landmark_range>                        Nearest_landmark_table_internal;

  typedef Landmark_id Vertex_handle;

 protected:
  Nearest_landmark_table_internal              nearest_landmark_table_;
  
 public:
  ActiveWitness* w_ref;
  bool first_spotted;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */

  //@{

  Witness_complex_cof()
  {
  }
  
  /**
   *  \brief Initializes member variables before constructing simplicial complex.
   *  \details Records nearest landmark table.
   *  @param[in] nearest_landmark_table needs to be a range of a range of pairs of nearest landmarks and distances.
 *         The range of pairs must admit a member type 'iterator'. The dereference type 
 *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'.
   */

  Witness_complex_cof(Nearest_landmark_table_ const & nearest_landmark_table)
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
                      std::size_t limit_dimension = std::numeric_limits<std::size_t>::max())    
  {
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;

    typedef std::map<Simplex_key, bool> Simplex_bool_map;
    
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

    Simplex_bool_map* prev_dim_map = new Simplex_bool_map();
    fill_vertices(max_alpha_square, complex, active_witnesses, prev_dim_map);

    Simplex_bool_map* dim1_map = new Simplex_bool_map();
    fill_edges(max_alpha_square, complex, active_witnesses, prev_dim_map, dim1_map);
    delete prev_dim_map;
    prev_dim_map = dim1_map;

    Landmark_id k = 2; /* current dimension in iterative construction */
    while (!active_witnesses.empty() && k <= limit_dimension) {
      typename ActiveWitness::Table::iterator it = w_ref->nearest_landmark_table_.begin();
      // std::advance(it, 2);
      // if (it->first != 3) {
      //   std::cout << "Attention!\n";
      // }

      Simplex_bool_map* curr_dim_map = new Simplex_bool_map();
      fill_simplices(max_alpha_square, k, complex, active_witnesses, prev_dim_map, curr_dim_map);

      delete prev_dim_map;
      prev_dim_map = curr_dim_map;
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
             typename SimplexBoolMap >
  void fill_vertices(const double alpha2,
                     SimplicialComplexForWitness& complex,
                     ActiveWitnessList& aw_list,
                     SimplexBoolMap* sw_map) const
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;

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
        complex.insert_simplex(Vertex_vector(1, l_it->first), filtration_value);
      }
    }
    for (auto aw_it = aw_list.begin(); aw_it != aw_list.end(); ++aw_it) {
      typename ActiveWitness::iterator l_it = aw_it->begin();
      typename ActiveWitness::iterator end = aw_it->end();
      double norelax_dist2 = std::numeric_limits<double>::infinity();
      for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
        Simplex_handle sh = complex.find(Vertex_vector(1, l_it->first));
        Siblings* sib = complex.self_siblings(sh);
        Vertex_handle v = sh->first;
        Simplex_key sk(sib,v);
        (*sw_map)[sk] = true; 
        aw_it->increase();
        if (l_it->second < norelax_dist2)
          norelax_dist2 = l_it->second;
      }
      // std::cout << aw_it->counter() << " "; 
    }
    std::cout << "0-dim active witness list size = " << aw_list.size() << "\n";
    // std::cout << "\n\n";
  }

  /* \brief Fills the map "edges -> witnesses for edges"
   */
  template < typename SimplicialComplexForWitness,
             typename ActiveWitnessList,
             typename SimplexWitnessMap >
  void fill_edges(double alpha2,
                  SimplicialComplexForWitness& complex,
                  ActiveWitnessList& aw_list,
                  SimplexWitnessMap* dim0_map,
                  SimplexWitnessMap* dim1_map)
  {
    Vertex_vector simplex;
    simplex.reserve(2);
    auto aw_it = aw_list.begin();
    // for (auto w: aw_list) {
    //   int num_simplices = 0;
    //   add_all_faces_of_dimension(1,
    //                              alpha2,
    //                              std::numeric_limits<double>::infinity(),
    //                              w.begin(),
    //                              simplex,
    //                              complex,
    //                              w.end(),
    //                              num_simplices,
    //                              dim1_map);
    //   assert(simplex.empty());
    // }
    while (aw_it != aw_list.end()) {
      int num_simplices = 0;
      add_all_faces_of_dimension(1,
                                 alpha2,
                                 std::numeric_limits<double>::infinity(),
                                 aw_it->begin(),
                                 simplex,
                                 complex,
                                 aw_it->end(),
                                 num_simplices,
                                 dim1_map);
      assert(simplex.empty());
      aw_it++;
    }

  // for (auto aw: aw_list)
    //   std::cout << aw.counter() << " "; 
    std::cout << "1-dim active witness list size = " << aw_list.size() << "\n";
  }

  /* \brief Fills the map "k-simplex -> witnesses for edges"
   */
  template < typename SimplicialComplexForWitness,
             typename ActiveWitnessList,
             typename SimplexBoolMap >
  void fill_simplices(double alpha2,
                      std::size_t k,
                      SimplicialComplexForWitness& complex,
                      ActiveWitnessList& aw_list,
                      SimplexBoolMap* prev_dim_map,
                      SimplexBoolMap* curr_dim_map)
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
    typedef Gudhi::Simplex_tree_vertex_subtree_iterator<SimplicialComplexForWitness> Vertex_subtree_iterator;
    typedef boost::iterator_range<Vertex_subtree_iterator> Vertex_subtree_range;

    // coface precalculation
    for (auto sw_pair: *prev_dim_map) {
      auto v_it = complex.simplex_vertex_range(sw_pair.first.simplex_handle()).begin();
      unsigned counter = 0; /* I need the first vertex before last */ 
      while (counter != k-2) {
        v_it++;
        counter++;
      }
      Vertex_handle v1 = *(v_it++);
      Vertex_handle v0 = *v_it;
      Vertex_subtree_range v0s_range(Vertex_subtree_iterator(&complex, v0, k-1),
                                    Vertex_subtree_iterator());
      Vertex_subtree_range v1s_range(Vertex_subtree_iterator(&complex, v1, k-1),
                                    Vertex_subtree_iterator());
      for (auto sh2: v0s_range) {
        Vertex_vector coface;
        if (check_if_neighbors(complex, sw_pair.first.simplex_handle(), sh2, coface)) {
          double filtration_value = 0;
          if (all_faces_in(coface, &filtration_value, complex)) {
            std::pair<Simplex_handle, bool> sh_bool = complex.insert_simplex(coface);
            if (sh_bool.second) {
              Siblings* sib = complex.self_siblings(sh_bool.first);
              Vertex_handle v = sh_bool.first->first;
              curr_dim_map->emplace(Simplex_key(sib,v), false);
            }
          }
        }
      }
      for (auto sh2: v1s_range) {
        Vertex_vector coface;
        if (check_if_neighbors(complex, sw_pair.first.simplex_handle(), sh2, coface)) {
          double filtration_value = std::numeric_limits<double>::infinity();
          if (all_faces_in(coface, &filtration_value, complex)) {
            std::pair<Simplex_handle, bool> sh_bool = complex.insert_simplex(coface);
            if (sh_bool.second) {
              Siblings* sib = complex.self_siblings(sh_bool.first);
              Vertex_handle v = sh_bool.first->first;
              curr_dim_map->emplace(Simplex_key(sib,v), false);
            }
          }
        }
      }
    }
    std::cout << "Candidates: " << curr_dim_map->size() << std::endl;
    auto aw_it = aw_list.begin();
    std::vector<Landmark_id> simplex;
    simplex.reserve(k+1);
    while (aw_it != aw_list.end()) {
      int num_simplices = 0;
      bool ok = add_all_faces_of_dimension(k,
                                           alpha2,
                                           std::numeric_limits<double>::infinity(),
                                           aw_it->begin(),
                                           simplex,
                                           complex,
                                           aw_it->end(),
                                           num_simplices,
                                           curr_dim_map);
      // std::cout << num_simplices << " ";
      assert(simplex.empty());
      if (!ok)
        aw_list.erase(aw_it++); //First increase the iterator and then erase the previous element
      else
        aw_it++;
    }
    // std::cout << "\n";
    std::cout << k << "-dim active witness list size = " << aw_list.size() << "\n";
    std::cout << "Candidates after: " << curr_dim_map->size() << std::endl;
    remove_non_witnessed_simplices(complex, curr_dim_map);
  }

  /* \brief Adds recursively all the faces of a certain dimension dim witnessed by the same witness.
   * Iterator is needed to know until how far we can take landmarks to form simplexes.
   * simplex is the prefix of the simplexes to insert.
   * The output value indicates if the witness rests active or not.
   */
  template < typename SimplicialComplexForWitness,
             typename SimplexBoolMap>
  bool add_all_faces_of_dimension(int dim,
                                  double alpha2,
                                  double norelax_dist2,
                                  typename ActiveWitness::iterator curr_l,
                                  std::vector<Landmark_id>& simplex,
                                  SimplicialComplexForWitness& complex,
                                  typename ActiveWitness::iterator end,
                                  int& num_simplices,
                                  SimplexBoolMap* curr_dim_map)
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
    typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForWitness::Siblings Siblings;
    typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;

    if (curr_l == end)
      return false;
    bool will_be_active = false;
    typename ActiveWitness::iterator l_it = curr_l;
    if (dim > 0)
      for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
        simplex.push_back(l_it->first);
        if (complex.find(simplex) != complex.null_simplex()) {
          typename ActiveWitness::iterator next_it = l_it;
          will_be_active = add_all_faces_of_dimension(dim-1,
                                                      alpha2,
                                                      norelax_dist2,
                                                      ++next_it,
                                                      simplex,
                                                      complex,
                                                      end,
                                                      num_simplices,
                                                      curr_dim_map) || will_be_active;
        }
        assert(!simplex.empty());
        simplex.pop_back();
        // If norelax_dist is infinity, change to first omitted distance
        if (l_it->second <= norelax_dist2)
          norelax_dist2 = l_it->second;
      } 
    else if (dim == 0)
      for (; l_it != end && l_it->second - alpha2 <= norelax_dist2; ++l_it) {
        simplex.push_back(l_it->first);
        double filtration_value = 0;
        // if norelax_dist is infinite, relaxation is 0.
        if (l_it->second > norelax_dist2) 
          filtration_value = l_it->second - norelax_dist2;
        if ((simplex[0] == 5) &&
            (simplex[1] == 7) && first_spotted) {
          w_ref = l_it.aw_;
          first_spotted = false;
        }
        // Two different modes: for edges and others
        if (simplex.size() == 2) {
          if (all_faces_in(simplex, &filtration_value, complex)) {
            will_be_active = true;
            auto sh_bool = complex.insert_simplex(simplex, filtration_value);
            if (sh_bool.second) {
              Simplex_handle sh = sh_bool.first;
              Vertex_handle v = sh->first;
              Siblings* sib = complex.self_siblings(sh);
              (*curr_dim_map)[Simplex_key(sib, v)] = true;
            }
            // if ((l_it.aw_->nearest_landmark_table_.begin()->first == 5) &&
            //     ((l_it.aw_->nearest_landmark_table_.begin()++)->first == 7) &&
            //     ((l_it.aw_->nearest_landmark_table_.begin()++++)->first == 3))
            //   w_ref = l_it.aw_;
            num_simplices++;
          }
        }
        // else {
        //   if (all_faces_in(simplex, &filtration_value, complex)) {
        //     will_be_active = true;
        //     auto sh_bool = complex.insert_simplex(simplex, filtration_value);
        //     if (sh_bool.second) {
        //       Simplex_handle sh = sh_bool.first;
        //       Vertex_handle v = sh->first;
        //       Siblings* sib = complex.self_siblings(sh);
        //       (*curr_dim_map)[Simplex_key(sib, v)] = true;
        //     }
        //     num_simplices++;
        //   }
        // }
        else {
          if ((simplex[0] == 5) &&
              (simplex[1] == 7) &&
              (simplex[2] == 9)) {
            w_ref = l_it.aw_;
            first_spotted = false;
          }
          Simplex_handle sh = complex.find(simplex);
          if (sh != complex.null_simplex()) {
            will_be_active = true;
            complex.insert_simplex(simplex, filtration_value);
            Vertex_handle v = sh->first;
            Siblings* sib = complex.self_siblings(sh);
            (*curr_dim_map)[Simplex_key(sib, v)] = true;
            num_simplices++;
          }
        }
        assert(!simplex.empty());
        simplex.pop_back();
        // If norelax_dist is infinity, change to first omitted distance
        if (l_it->second < norelax_dist2)
          norelax_dist2 = l_it->second;
      }
    return will_be_active;
  }

  /** Removes the non-witnessed candidate simplices.
   */
  template < class SimplicialComplexForWitness,
             class SimplexBoolMap >
  void remove_non_witnessed_simplices(SimplicialComplexForWitness& complex, SimplexBoolMap* curr_dim_map)
  {
    std::list<typename SimplexBoolMap::key_type> elements_to_remove;
      
    for (auto sw: *curr_dim_map) {
      if (!sw.second) {
        complex.remove_maximal_simplex(sw.first.simplex_handle());
        elements_to_remove.push_back(sw.first);
      }        
    }
    for (auto key: elements_to_remove)
      curr_dim_map->erase(key);
  }

  template < class SimplicialComplexForWitness>
  void print_complex(SimplicialComplexForWitness& complex)
  {
    std::cout << complex << std::endl;
  }
  
};
  
}  // namespace witness_complex

}  // namespace Gudhi

#endif  // WITNESS_COMPLEX_COF_H_
