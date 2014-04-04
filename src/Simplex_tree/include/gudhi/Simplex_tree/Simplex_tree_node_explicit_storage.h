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

/**
 * \brief Node of a simplex tree with filtration value
 * and simplex key.
 *
 * It stores explicitely its own filtration value and its own Simplex_key.
 */
template < class SimplexTree > 
class Simplex_tree_node_explicit_storage {
  public:
  friend SimplexTree;

  typedef typename SimplexTree::Siblings         Siblings;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Simplex_key      Simplex_key;

  //private:
  //friend class Simplex_tree; 
  // Default constructor.
  Simplex_tree_node_explicit_storage() :
  children_(NULL),
  simplex_key_(-1),
  filtration_(0) {}

 /* Simplex_tree_node_explicit_storage( Filtration_value filtration ) :
  children_(NULL),
  simplex_data_(0),
  filtration_(filtration) 
  {}*/


  Simplex_tree_node_explicit_storage(Siblings * sib,
                                     Filtration_value filtration) :
  children_(sib),
  simplex_key_(-1),
  filtration_(filtration) {}


  void assign_key(Simplex_key key) { simplex_key_ = key; }

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
  void assign_children (Siblings * children) { children_ = children; }
  /**
   *
   */
  void assign_filtration(double filtration_value) { filtration_ = filtration_value; }
  
  Filtration_value filtration() { return filtration_; }

  /** Careful -> has_children() must be true*/
  Siblings * children() { return children_; }
  
  Simplex_key key() { return simplex_key_; }
  
private:  
   Siblings *              children_;
  
  // Data attached to simplex, explicit storage
  Simplex_key             simplex_key_;
  Filtration_value        filtration_;   //value in the filtration
  
};

#endif // GUDHI_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H
