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
         , class N
         , class MC >
class Simplex_tree_siblings;

/**
 * \brief Node of a simplex tree with filtration value
 * and simplex data.
 *
 * It stores explicitely its own filtration value and its own Simplex_data.
 */
template < class SimplexTree > 
class Simplex_tree_node_explicit_storage {
  public:
friend SimplexTree;

  typedef typename SimplexTree::Siblings         Siblings;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Simplex_data     Simplex_data;

  //private:
  //friend class Simplex_tree; 
  // Default constructor.
  Simplex_tree_node_explicit_storage() :
  children_(NULL),
  simplex_data_(-1),
  filtration_(0) {}

 /* Simplex_tree_node_explicit_storage( Filtration_value filtration ) :
  children_(NULL),
  simplex_data_(0),
  filtration_(filtration) 
  {}*/


  Simplex_tree_node_explicit_storage(Siblings * sib,
                                     Filtration_value filtration) :
  children_(sib),
  simplex_data_(-1),
  filtration_(filtration) {}


  void assign_data(Simplex_data key) { simplex_data_ = key; }

  /**
   * Return true if the node has children,
   * false otherwise.
   */
  //bool has_children(Vertex label)
  //{ //if(children_ == NULL)             return false; //for root simplices
  //  return (children_->parent() == label);}
  /**
   * Assign a children to the node
   */
  void 
  assign_children (Siblings *      children)
  { children_ = children; }
  /**
   *
   */
  void assign_filtration(double filtration_value)
  {  filtration_ = filtration_value;  }
  
  Filtration_value filtration()
  { return filtration_; }

  /** Careful -> has_children() must be true*/
  Siblings *     children()
  { return children_; }
  
  Simplex_data simplex_data()
  { return simplex_data_; }
  
  Simplex_data data() { return simplex_data_; }

private:  
   Siblings *              children_;
  
  // Data attached to simplex, explicit storage
  Simplex_data             simplex_data_;
  Filtration_value         filtration_;   //value in the filtration
  
};

#endif // GUDHI_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H
