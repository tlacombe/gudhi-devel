/*
 *  Simplex_tree_node_explicit_storage.cpp
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#include "Simplex_tree_siblings.h"

/*template <class V, class F, class SD>
Simplex_tree_node_explicit_storage<V,F,SD>::
Simplex_tree_node_explicit_storage() :
  children_(NULL),simplex_data_(0),filtration_(0) {}
*/
/*template <class V, class F, class SD>
Simplex_tree_node_explicit_storage<V,F,SD>::  
Simplex_tree_node_explicit_storage( F filtration,
                            			  Simplex_tree_siblings< V,
                                                  Simplex_tree_node_explicit_storage//<V,F,SD>
                                                         >     *self_sib) :
  children_(self_sib),
  simplex_data_(0),
  filtration_(filtration) 
{}
*/
/*
template <class V, class F, class SD>
Simplex_tree_node_explicit_storage<V,F,SD>::  
Simplex_tree_node_explicit_storage(F filtration) :
  children_(NULL),
  simplex_data_(0),
  filtration_(filtration) 
{}*/
/*  
template <class V, class F, class SD>
Simplex_tree_node_explicit_storage<V,F,SD>::  
Simplex_tree_node_explicit_storage(const Simplex_tree_node_explicit_storage &other) :
  children_(other.children_),
  simplex_data_(other.simplex_data_),
  filtration_(other.filtration_)  
{}*/

/** \todo remove next_sib == NULL test in Node.cpp*/
template <class V, class F, class SD>
Simplex_tree_siblings< V,
                       F,
                       Simplex_tree_node_explicit_storage<V,F,SD> 
                     > *
Simplex_tree_node_explicit_storage<V,F,SD>::
self_siblings(V label)
{
  Simplex_tree_siblings< V,
                         F,
                         Simplex_tree_node_explicit_storage<V,F,SD>
                       > * next_sib = children_;
  if(next_sib == NULL) std::cerr << "Error in self_siblings \n";
  if(next_sib->parent() == label) return next_sib->oncles();
  else                            return next_sib;
}

/*template <class V, class F, class SD>
void
Simplex_tree_node_explicit_storage<V,F,SD>::
list_of_vertices( std::vector< V > &sigma,
						      V             label)
{
  Simplex_tree_siblings< Filtered_simplex_tree_node<S,V,F> > * curr_sib =
    self_siblings(label);
  sigma.push_back(label);
  
  while(curr_sib != NULL)
    {
      sigma.push_back(curr_sib->parent());
      curr_sib = curr_sib->oncles();
    }
}*/

template <class V, class F, class SD>
bool
Simplex_tree_node_explicit_storage<V,F,SD>::
is_edge(V label)
{
  if(self_siblings(label)->oncles() == NULL) return true;
  else return false;
}

/*template <class V, class F, class SD>
bool 
Simplex_tree_node_explicit_storage<V,F,SD>::
has_children(V label)
{
  if(children_ == NULL)             return false; //for root simplices
  if(children_->parent() == label)  return true;
  else                              return false;
}*/

template <class V, class F, class SD>
int 
Simplex_tree_node_explicit_storage<V,F,SD>::
dimension(V label)
{
  Simplex_tree_siblings< V,
                         F,
                       Simplex_tree_node_explicit_storage<V,F,SD>
                       > * curr_sib = self_siblings(label);
  int dim = 0;
  while(curr_sib != NULL)
    { ++dim;
      curr_sib = curr_sib->oncles(); }
  return dim;
}
