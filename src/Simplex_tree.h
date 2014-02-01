/*
 *  Simplex_tree.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_SIMPLEX_TREE_H
#define GUDHI_SIMPLEX_TREE_H

#include <boost/container/flat_map.hpp>

#include "Simplex_tree_node_explicit_storage.h"
#include "Simplex_tree_siblings.h"

#include "Simplex_tree_iterators.h" // implementation of the iterators
                                    // for Simplex_tree.

/**
 * \brief Simplex tree data structure.
 *
 * \implements SimplexDataFilteredSimplicialComplexDS
 */
template < class MetricSpace >//, class SimplexTreeNode , class SimplexTreeSiblings >
class Simplex_tree {
  public:
  /// \name Type definitions
  /// @{
  typedef typename MetricSpace::FT                    Filtration_value     ;
  typedef typename MetricSpace::Vertex                Vertex               ;   //-
  typedef typename MetricSpace::Space_vertex_iterator Space_vertex_iterator;
  typedef typename MetricSpace::Space_vertex_range    Space_vertex_range   ;

  /** \brief Type of data store in each simplex. */
  typedef int                                         Simplex_data;
  /** \brief Node in the simplex tree.*/ 
  typedef Simplex_tree_node_explicit_storage < Vertex,
                                               Filtration_value,
                                               Simplex_data >    Node          ;
 /** \brief Must be an ordered range. */
 // typedef boost::flat_map< Vertex, Node >                        Dictionary    ;

  /** \brief Set of nodes sharing a same parent in the simplex tree. */
  typedef Simplex_tree_siblings < Vertex
                                , Filtration_value
                                , Node >            Siblings      ;
/** \todo Didn't manage to define the Dictionary type in Simplex_tree because
* of the mutual dependence of Siblings and Node.*/
 typedef typename Siblings::Dictionary              Dictionary;
 typedef typename Dictionary::iterator Dictionary_it;
/**
  * \todo typedef   Simplex_handle
  */
  typedef typename Dictionary::iterator                      Simplex_handle;


  /// @}

  /// \name Simplex Vertex traversal    
  /// @{
 /**
   * Returns a Simplex_vertex_range for the sequence of all
   * vertices of the simplex associated
   * to Simplex_handle sh.
   */
  Simplex_vertex_range< Simplex_tree > simplex_vertex_range(Simplex_handle sh)
  { return Simplex_vertex_range< Simplex_tree > (sh); }

  void print(Simplex_handle sh, std::ostream& os = std::cout)
  { Simplex_vertex_range< Simplex_tree > svr = simplex_vertex_range(sh);
    for(Simplex_vertex_iterator< Simplex_tree > it = svr.begin();
        it != svr.end(); ++it) {os << *it << " ";}
    os << std::endl;}
  /// @}


  


  /** \brief Empty constructor.*/
  Simplex_tree( MetricSpace & ms) :
  ms_(&ms),
  nb_vertices_(0),
  size_cpx_(0),
  root_(),
  filtration_vect_()
  //NULL_sh_(ms_->NULL_vertex,Node())
  {}


/// \name Simplex container methods.
/// @{

/** 
* \brief Given a sequence of Vertices, returns the
* Simplex_handle in the simplex tree corresponding 
* to the simplex with this set of Vertices.
* 
* The sequence of Vertices must be sorted in 
* increasing order.
*
* If the simplex is not in the simplex tree, returns end().
*
* \todo Simplex_tree find and insert.
*/

/*Simplex_handle find(std::vector < Vertex > & s)
{ if(s.size() == 0) std::cerr << "Empty simplex \n";
  if(s.size() == 1) std::cerr << "Vertex \n";

  typename std::vector< Vertex >::iterator it = s.begin();
  if(! root_[ *it ].has_children( *it )) { }    // not there
  Siblings * for_sib = root_[ *it ].children();
  ++it;
  Simplex_handle sh = for_sib->find( *it ); //some stop condition here
  ++it;
  for( ; it != s.end(); ++it)
    { if(! sh->second.has_children(sh->first)) std::cerr << "Not here \n";    //must create a st.end() Simplex_handle...
      sh = sh->second.children()->find(*it); }// return some false if not here...
  return sh; }
   */

//Simplex_handle insert(); //input a vertex_range


/**
 * Returns a Complex_simplex_range over all simplices
 * of dimension > 0 
 * in the simplicial complex stored in the simplex
 * tree.
 */
Complex_simplex_range< Simplex_tree > complex_simplex_range()
{ return Complex_simplex_range< Simplex_tree > (this); }
  /// @}

  
  /// \name Acces methods
  /// @{
  /** Returns a pointer to the geometry traits.*/
  MetricSpace *                     ms()                { return ms_; }
  
  /** Returns the maximal threshold for the Filtration_value.*/
 // Filtration_value                rho_max()           { return rho_max_; }
  /** Returns the number of vertices in the complex.*/
  size_t                            nb_vertices()       { return nb_vertices_; }
  /** \brief Returns the number of faces of the complex.*/
  size_t                            size_complex()      { return size_cpx_; }
  /** Returns a reference to the root nodes of the simplex tree.*/
  std::vector< Node > &             root()              { return root_; }
  /** Returns the vector representing the filtration order.*/
  //std::vector< Simplex_handle > & filtration_vector() { return filtration_vect_; }
  /// @}
       
  /// \name Boundary Simplex traversal
  /// @{
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
  Boundary_simplex_range< Simplex_tree > boundary_simplex_range(Simplex_handle sh)
  { return Boundary_simplex_range< Simplex_tree > (this,sh);}


/**
 * \brief Iterator over the simplices of a filtration.
 *
 * 'value_type' must be Simplex_handle.
 */
typedef typename std::vector< Simplex_handle >::iterator Filtration_simplex_iterator;  
/**
 * \brief Range over the simplices following the order
 * of a filtration.
 *
 * Methods .begin() and .end() return
 * a Filtration_simplex_iterator. 
 */
typedef typename std::vector< Simplex_handle > Filtration_simplex_range;
/**
* \brief Returns a Filtration_simplex_range for the sequence of all
* simplices of the simplicial complex, in the order of a filtration.
*/
Filtration_simplex_range filtration_simplex_range() 
{ if(filtration_vect_.size() == 0) { initialize_filtration(); }
  return Filtration_simplex_range(this);}
/**
* \brief Initializes the filtration.
*
* Will be automatically called when calling
* filtration_simplex_range() if the filtration has
* not been initialized yet.
*/ 
void initialize_filtration();
/**
* Returns true iff sh1 is a subface of sh2.
*/
bool is_subface(Simplex_handle sh1, Simplex_handle sh2);
/**
* The use of stable_sort + the is_subface comparison
* allows a nice filtration strategy order.
*/
struct is_before_in_filtration {
  is_before_in_filtration(Simplex_tree * st) :
  st_(st) {}

  bool operator()(const Simplex_handle sh1,
                  const Simplex_handle sh2) const 
  { if(sh1->second.filtration() != sh2->second.filtration())
    {return sh1->second.filtration() < sh2->second.filtration();}
    return !(st_->is_subface(sh1,sh2)); }  //is sh1 a subface of sh2
    
  Simplex_tree * st_;
};

/**
 * \brief Inserts a 1-skeleton.
 *
 * Inserts all edges given by Neighbor_vertex_range
 * defined in NeighborsGeometryTraits, which produces, for a given
 * vertex, a range allowing to iterate over its neighbors in 
 * the 1-skeleton.
 *
 * \todo Works on a range of points. + easily parallelizable.
 * attention static vector => not thread safe in intersection.
 */
template< class NeighborGraph >
void insert_graph( NeighborGraph & ng )
{
  root_ = std::vector< Node >( ng.size_graph() , Node() );
  Space_vertex_range v_rg = ms_->space_vertex_range();
  for(Space_vertex_iterator v_it = v_rg.begin();
      v_it != v_rg.end(); ++v_it)
  {
    typename NeighborGraph::adjacency_range n_range =
                                   ng.vertex_adjacency_range(*v_it);
    for(typename NeighborGraph::adjacency_iterator n_it = n_range.begin();
        n_it != n_range.end(); ++n_it)
    {
      if(*v_it < *n_it)       //count edges only once. 
        {
          if(! root_[*v_it].has_children(*v_it)) 
          { this->root_[*v_it].assign_children(new Siblings(NULL,*v_it)); }
          this->root_[*v_it].children()->insert(*n_it,
                                                ms_->distance(*v_it,*n_it)
                                                );
        }
    }
  }
// Update size of the complex
this->size_cpx_ += this->root_.size();
int v = 0;
for(typename std::vector< Node >::iterator r_it = root_.begin();
    r_it != root_.end(); ++r_it,++v)
  { if(r_it->has_children(v)) { size_cpx_ += r_it->children()->members().size();}}
}

/**
  * \brief Expansion of the simplicial complex until
  * dimension max_dim.
  *
  * \todo Define expansion formally; what if more than a graph 
  * in the complex?
  */
void expansion(int max_dim)
{
  int curr_vertex = 0;
  for(typename std::vector< Node >::iterator root_it = root_.begin();
      root_it != root_.end(); ++root_it, ++curr_vertex)
  {
    if(root_it->has_children(curr_vertex)) 
    { siblings_expansion(root_it->children(), max_dim-1); }
  }
}
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
                  Filtration_value                             filtration)
{
  if(begin1 == end1 || begin2 == end2) return;// 0;
  while( true ) 
  {
    if( begin1->first == begin2->first )
    {
      intersection.push_back(std::pair< Vertex, Node >(begin1->first,
                             Node(NULL,
                                  maximum(begin1->second.filtration(),
                                          begin2->second.filtration(),
                                          filtration))));
      ++begin1;  ++begin2;
      if( begin1 == end1 || begin2 == end2 ) return;
    }
    else 
    { 
      if( begin1->first < begin2->first ) 
      { ++begin1;
        if(begin1 == end1) return; }
      else 
      { ++begin2;
        if(begin2 == end2) return;
      }
    }
  }
}

private:    
  /** Maximum over 3 values.*/
  Filtration_value maximum(Filtration_value a, 
                           Filtration_value b, 
                           Filtration_value c )
  { Filtration_value max = ( a < b ) ? b : a;
    return ( ( max < c ) ? c : max ); }

/** \brief Pointer to a metric space. */
  MetricSpace *                    ms_             ;
/** \brief Threshold for the filtration function. */
//  Filtration_value                 rho_max_        ;
/** \brief Number of vertices. The set of vertices is static.*/      
  int                              nb_vertices_    ;
/** \brief Total number of simplices in the complex, without the empty simplex.*/
  int                              size_cpx_       ;
/** \brief Set of simplex tree Nodes representing the vertices.*/  
  std::vector< Node >              root_           ;
/** \brief Simplices ordered according to a filtration*/  
  std::vector< Simplex_handle >    filtration_vect_;
/** \brief A NULL Simplex_handle; useful for the implementation.*/  
//    Simplex_handle                   NULL_sh_        ;

};



#include "Simplex_tree.hpp"         //implementation of the methods in
                                    // Simplex_tree.

#endif // GUDHI_FLAG_SIMPLEX_TREE_H
