/*
 *  Simplex_tree_node_explicit_storage.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H
#define GUDHI_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H

#include <vector>
#include <iostream>

template < class V
         , class F
         , class N >
class Simplex_tree_siblings;

/**
 * \brief Node of a simplex tree with filtration value
 * and simplex data.
 *
 * It stores explicitely its own filtration value and its own Simplex_data.
 */
template < class Vertex
         , class Filtration_value
         , class Simplex_data = int > 
class Simplex_tree_node_explicit_storage {
  public: 
  /**
   * Default constructor.
   */
  Simplex_tree_node_explicit_storage() :
  children_(NULL),
  simplex_data_(-1),
  filtration_(0) {}

 /* Simplex_tree_node_explicit_storage( Filtration_value filtration ) :
  children_(NULL),
  simplex_data_(0),
  filtration_(filtration) 
  {}*/

//template < class V, class F, class N >
  Simplex_tree_node_explicit_storage( 
                Simplex_tree_siblings< Vertex
                                     , Filtration_value
                                     , Simplex_tree_node_explicit_storage > * sib,
                Filtration_value                                        filtration) :
  children_(sib),
  simplex_data_(-1),
  filtration_(filtration) {}

  /**
   * Constructor with values.
   */
  // Simplex_tree_node_explicit_storage( Filtration_value filtration,
  //                                     Simplex_tree_siblings< Vertex
  //                                                          , Simplex_tree_node_explicit_storage
  //                                                          >                *self_sib);
   /**
   * Constructor with values.
   */

  /** Necessary
   * Copy constructor
   */
  //Simplex_tree_node_explicit_storage(const Simplex_tree_node_explicit_storage &other);
  /**
   * Destructor
   */
  //~Simplex_tree_node_explicit_storage();
  /**
   * When in a simplex tree, returns a pointer
   * to the Simplex_tree_siblings containing the node.
   *
   * Vertex label is the biggest Vertex in the simplex
   * represented by the node.
   */
  //template <class N>
  Simplex_tree_siblings< Vertex
                       , Filtration_value
                       , Simplex_tree_node_explicit_storage > * 
  self_siblings(Vertex label);
  /**
   * Return true if the node has children,
   * false otherwise.
   */
  bool has_children(Vertex label)
  { if(children_ == NULL)             return false; //for root simplices
    if(children_->parent() == label)  return true;
    else                              return false;}
  /**
   * Assign a children to the node
   */
  void 
  assign_children( Simplex_tree_siblings< Vertex
                                        , Filtration_value
                                        , Simplex_tree_node_explicit_storage > *      children)
  { children_ = children; }
  /**
   *
   */
  void assign_filtration(double filtration_value)
  {  filtration_ = filtration_value;  }
  /**
   * Return the dimension of the simplex corresponding
   * to the node.
   */
  int dimension(Vertex label);
  /**
   * Fill sigma with the list of vertices
   * of the simplex corresponding to the node
   */    
  //void list_of_vertices( std::vector< Vertex > & sigma,
  //                       Vertex                  label );
  /**
   * Return true if the simplex corresponding 
   * to the node is an edge, false otherwise.
   */
  bool is_edge( Vertex label );
  
  Filtration_value filtration()
  { return filtration_; }

  /** Careful -> has_children() must be true*/
  Simplex_tree_siblings< Vertex
                       , Filtration_value
                       , Simplex_tree_node_explicit_storage > *     children()
  { return children_; }
  
  Simplex_data simplex_data()
  { return simplex_data_; }
  

  /***************************/
private:  
   Simplex_tree_siblings< Vertex
                       , Filtration_value
                       , Simplex_tree_node_explicit_storage > *     children_;
  //S_inter_Node        * inter_node_;
  //void            * inter_node_;
  
  // Data attached to simplex;
  Simplex_data             simplex_data_;
  Filtration_value         filtration_; //value in the filtration
  
};

#endif // GUDHI_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H
