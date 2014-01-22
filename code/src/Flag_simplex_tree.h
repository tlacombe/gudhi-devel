/*
 *  Rips_simplex_tree.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_FLAG_SIMPLEX_TREE_H
#define GUDHI_FLAG_SIMPLEX_TREE_H

#include <boost/container/flat_map.hpp>
#include "Euclidean_rips_naive_geometry_traits.h"
#include "Filtered_simplex_tree_node.h"
#include "Simplex_tree_siblings.h"

/**
 * \brief Simplex tree data structure to construct flag complexes
 *
 * The type of complex is contained in the template
 * parameter NeighborsGeometryTraits which furnishes 
 * a range for the neighbors in the graph of a vertex
 *
 * \implements SimplexDataFilteredSimplicialComplexDS
 */
template < class NeighborsGeometryTraits >
class Flag_simplex_tree {
 public:
  /// \name Type definitions
  /// @{
  typedef typename NeighborsGeometryTraits::FT              Filtration_value          ;
  typedef typename NeighborsGeometryTraits::Point           Point                     ;
  typedef typename NeighborsGeometryTraits::Point_range     Point_range               ;
  typedef typename NeighborsGeometryTraits::Point_iterator  Point_iterator            ;
  typedef typename NeighborsGeometryTraits::Vertex          Vertex                    ;   //-
  typedef typename NeighborsGeometryTraits::Vertex_iterator Vertex_iterator           ;
  typedef typename NeighborsGeometryTraits::Vertex_range    Vertex_range              ;
  typedef typename NeighborsGeometryTraits::Neighbor_vertex_range 
                                                            Neighbor_vertex_range     ;
  typedef typename NeighborsGeometryTraits::Neighbor_vertex_iterator 
                                                            Neighbor_vertex_iterator  ;
  /** Node in the simplex tree.*/ 
  typedef          Filtered_simplex_tree_node               Node                      ;
  /** \brief Set of nodes sharing a same parent in the simplex tree. */
  typedef          Simplex_tree_siblings                    Siblings                  ;
  /** \brief Must be an ordered range. */
  typedef typename Siblings::Dictionary                     Dictionary                ;
  /** \todo Probably not correct, should be Dictionary::iterator. */
  typedef typename Siblings::Dictionary_it                  Dictionary_it             ;

  typedef typename Siblings::Dictionary_it                  Simplex_handle            ;
  /// @}


  /// \name Simplex Vertex traversal
  /// @{
  /**
   * \brief Range over the vertices of a simplex
   *
   * Methods .begin() and .end() return
   * a Simplex_vertex_iterator. 
   *
   * The order in which the Vertices are visited defines
   * the canonical orientation of the simplex.
   */
  class              Simplex_vertex_range        ;
  /**
   * \brief Iterator over the vertices of a simplex
   *
   * 'value_type' must be Vertex.
   */
  class              Simplex_vertex_iterator     ;
  /**
   * Returns a Simplex_vertex_range for the sequence of all
   * vertices of the simplex associated
   * to Simplex_handle sh.
   */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle sh)
  { return Simplex_vertex_range(sh); }

  void print(Simplex_handle sh, std::ostream& os = std::cout);
  /// @}


  /** \brief Empty constructor.*/
  Flag_simplex_tree();

  /** Destructor.*/
  ~Flag_simplex_tree()  { delete gt_; }


  /// \name Simplex container methods.
  /// @{
  /** \todo svr must be ordered */
  /*    Simplex_handle does_simplex_belong_to_complex(std::vector< Vertex > & svr)
    {    
    if(svr.begin() == svr.end()) std::cerr << "Empty simplex \n";

    Simplex_vertex_iterator it = svr.begin();
    Simplex_handle sh = root_[*it];
    ++it;
    for(;    it != svr.end(); ++it)
    {
    if(! sh->second.has_children(sh->first))
    return std::cerr << "Not here \n";       //must create a st.end() Simplex_handle...
    sh = sh->second.children()->find(*it);            // return some false if not here...
    }
    return sh;
    }
  */

/** 
* \brief Given a sequence of Vertices, returns the
* Simplex_handle in the simplex tree corresponding 
* to the simplex with this set of Vertices.
* 
* The sequence of Vertices must be sorted in 
* increasing order.
*
* If the simplex is not in the simplex tree, returns end().
*/
Simplex_handle find(std::vector < Vertex > & s);

   /**
   * \brief Range over the simplices of a simplicial complex.
   *
   * Methods .begin() and .end() return
   * a Complex_simplex_iterator. 
   */
  class Complex_simplex_range;
  /**
   * \brief Iterator over the simplices of a 
   * simplicial complex.
   *
   * 'value_type' must be Simplex_handle.
   */
  class Complex_simplex_iterator;
  /**
   * Returns a Complex_simplex_range over all simplices
   * of dimension > 0 
   * in the simplicial complex stored in the simplex
   * tree.
   */
  Complex_simplex_range complexe_simplex_range()
  {return Complex_simplex_range(this);}
  /// @}

  /// \name Construction of the flag complex
  /// @{
  /**
   * \brief Construct the flag complex.
   *
   * First inserts all edges given by Neighbor_vertex_range
   * defined in NeighborsGeometryTraits, which produces, for a given
   * vertex, a range allowing to iterate over its neighbors in 
   * the 1-skeleton.
   * Then, realizes the expansion of the graph until dimension
   * dim_max.
   *
   * \todo Works on a range of points.
   */
   //       Point_range_sc  & point_range,
  void init(int              dim_max,
            Filtration_value rho_max);
  /**
   * Recursive expansion of the simplex tree.
   */
  void siblings_expansion(Siblings * siblings, //must contain elements
                          int k);
  /**
   * Intersects Dictionary 1 [begin1;end1) with
   * Dictionary 2 [begin2,end2) and assigns the
   * maximal possible Filtration_value to the Nodes.
   */
  void intersection(std::vector< std::pair< Vertex, Node > > &   intersection,
                    Dictionary_it                                begin1,
                    Dictionary_it                                end1,
                    Dictionary_it                                begin2,
                    Dictionary_it                                end2,
                    Filtration_value                             filtration);
  /// @}

  /// \name Acces methods
  /// @{
  /** Returns a pointer to the geometry traits.*/
  NeighborsGeometryTraits      *  gt()                { return gt_; }
  /** Returns the maximal threshold for the Filtration_value.*/
  Filtration_value                rho_max()           { return rho_max_; }
  /** Returns the number of vertices in the complex.*/
  int                             nb_vertices()       { return nb_vertices_; }
  /** \brief Returns the number of faces of the complex.*/
  int                             size_complex()      { return size_cpx_; }
  /** Returns a reference to the root nodes of the simplex tree.*/
  std::vector< Node > &           root()              { return root_; }
  /** Returns the vector representing the filtration order.*/
  std::vector< Simplex_handle > & filtration_vector() {return filtration_vect_; }

  /// @}
       
  /// \name Boundary Simplex traversal
  /// @{
   /**
   * \brief Range over the simplices in the boundary of a simplex.
   *
   * Methods .begin() and .end() return a Boundary_simplex_iterator.
   */
  class Boundary_simplex_range;
  /**
   * \brief Iterator over the simplices of the boundary of a
   * simplex.
   *
   * `value_type` must be `Simplex_handle`.
   */
  class Boundary_simplex_iterator;
  /**
   * \brief Returns a range over all simplices of the boundary of a simplex, i.e.
   * the set of codimension $1$ subsimplices of the simplex.
   *
   * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
   * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
   * simplices of the boundary in the order: 
   * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
   *
   * We note that the alternate sum of the simplices given by the iterator
   * gives the chains corresponding to the boundary of the simplex.
   */
  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh)
  { return Boundary_simplex_range(this,sh);}
  /// @}

  /// \name Filtration Simplex traversal
  /// @{
 /**
   * \brief Range over the simplices following the order
   * of a filtration.
   *
   * Methods .begin() and .end() return
   * a Filtration_simplex_iterator. 
   */
  class Filtration_simplex_range;
/**
   * \brief Iterator over the simplices of a filtration.
   *
   * 'value_type' must be Simplex_handle.
   */
  typedef std::vector< Simplex_handle >::iterator Filtration_simplex_iterator;
/**
   * \brief Returns a Filtration_simplex_range for the sequence of all
   * simplices of the simplicial complex, in the order of a filtration.
   */
  Filtration_simplex_range filtration_simplex_range() 
  {
    if(filtration_vect_.size() == 0) {initialize_filtration();}
    return Filtration_simplex_range(this);
  }
/**
  * \brief Initializes the filtration.
  *
  * Will be automatically called when calling
  * filtration_simplex_range() if the filtration has
  * not been initialized yet.
  */ 
  void initialize_filtration();
  /**
   * The use of stable_sort + the is_subface comparison
   * allows a nice filtration strategy order.
   */
  bool compare_simplices_fil(const Simplex_handle sh1,
                             const Simplex_handle sh2);
  /**
   * Returns true iff sh1 is a subface of sh2.
   */
  bool is_subface(Simplex_handle sh1, Simplex_handle sh2);
  /// @}
  

 private:    
  /**
   * Maximum over 3 values.
   */
  Filtration_value maximum(Filtration_value a, 
               Filtration_value b, 
               Filtration_value c )
  {
    Filtration_value max = ( a < b ) ? b : a;
    return ( ( max < c ) ? c : max );
  }


  NeighborsGeometryTraits        *    gt_            ;
  Filtration_value                rho_max_        ;    
  int                            nb_vertices_;
  int                            size_cpx_        ; //with or without vertices ?
  //    int                            dimension_cpx_    ; 
  std::vector< Node >                root_            ;  //set of top nodes

  std::vector< Simplex_handle >     filtration_vect_;

  Simplex_handle                    NULL_sh_        ;
};










std::ostream& operator<<(std::ostream& os, Simplex_tree_siblings & obj)
{
  os << "--Oncles: @ " << (long int)(obj.oncles()) << "\n";
  os << "--Parent:   " << obj.parent() << "\n";
  os << "Siblings: @ " << (long int)(&obj) << "\n";
  for(Simplex_tree_siblings::Dictionary_it sh = obj.members().begin();
      sh != obj.members().end(); ++sh)
    {    os << "[" << sh->first << ":" << sh->second.filtration() <<"] ";    }
  os << std::endl << std::endl;
  for(Simplex_tree_siblings::Dictionary_it sh = obj.members().begin();
      sh != obj.members().end(); ++sh)
    {if(sh->second.has_children(sh->first)) os << *(sh->second.children());}
  return os;
};
std::ostream& operator<<(std::ostream& os, Flag_simplex_tree< Euclidean_rips_naive_geometry_traits > & obj)
{
  os << "Flag Simplex Tree: \n";
  os << "Size Cpx   = " << obj.size_complex() << std::endl;
  //    os << "Dimension  = " << obj.dimension_cpx_ << std::endl;
  os << "rho_max = " << obj.rho_max() << std::endl;
  os << "nb_V    = " << obj.nb_vertices() << std::endl;
  os << std::endl;

  int v = 0;
  os << "@ 0000000000:   ";
  for(std::vector< Filtered_simplex_tree_node >::iterator it = obj.root().begin();
      it != obj.root().end(); ++it, ++v)
    {    os << v << " ";}
  os << std::endl;

  v = 0;
  for(std::vector< Filtered_simplex_tree_node >::iterator it = obj.root().begin();
      it != obj.root().end(); ++it,++v)
    {
      if(it->has_children(v)) os << *(it->children());    
    }
  return os;
};    


#include "Flag_simplex_tree_iterators.h" 
#include "Flag_simplex_tree.hpp"


#endif // GUDHI_FLAG_SIMPLEX_TREE_H
