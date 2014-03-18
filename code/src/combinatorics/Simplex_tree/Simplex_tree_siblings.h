/*
 *  Simplex_tree_siblings.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_SIMPLEX_TREE_SIBLINGS
#define GUDHI_SIMPLEX_TREE_SIBLINGS

#include "boost/container/flat_map.hpp"
#include "Simplex_tree_node_explicit_storage.h"

/**
 * 
 */
template < 
// class Vertex
//          , class Filtration_value
//          , class Node
          class SimplexTree
         , class MapContainer >// < Vertex, Node >
       //   >//, class Dictionary > 
class Simplex_tree_siblings {
public:
  friend SimplexTree;

  typedef typename SimplexTree::Vertex Vertex;
  typedef typename SimplexTree::Filtration_value  Filtration_value;
  typedef typename SimplexTree::Node Node;
  typedef  MapContainer                           Dictionary;
  typedef typename MapContainer::iterator         Dictionary_it;

//typedef boost::container::flat_map< Vertex, Node > Dictionary;
  
  // Default constructor
  Simplex_tree_siblings() 
  : oncles_(NULL)
  , parent_(-1)
  , members_() {}
  
  // Construct with values
  Simplex_tree_siblings(Simplex_tree_siblings  * oncles,
                        Vertex                   parent )
  : oncles_(oncles)
  , parent_(parent)
  , members_() {}
  
  /**
   * 'members' is sorted and unique.
   */
  Simplex_tree_siblings(Simplex_tree_siblings * oncles,
                        Vertex                  parent,
                        std::vector< std::pair< Vertex, Node > > & members) :
  oncles_(oncles),
  parent_(parent),
  members_(boost::container::ordered_unique_range,members.begin(),members.end())
  {
    for(auto map_it = members_.begin();
        map_it != members_.end(); map_it++)
    {  map_it->second.assign_children(this);  }
  }
  
  

  /**
   * Construct with initialized set of members
   */
  /*Simplex_tree_siblings(Simplex_tree_siblings * oncles,
                        Vertex                  parent, 
                        Dictionary            & init_members) :
  oncles_(oncles),
  parent_(parent)
  {
    members_ = Dictionary(init_members);
    for(auto map_it = members_.begin();
        map_it != members_.end(); map_it++)
    {  map_it->second.assign_children(this);  }
  }*/
    
  
  /**
   * \brief Inserts a Node in the set of siblings nodes.
   *
   * If already present, assigns the minimal filtration value 
   * between input filtration_value and the value already 
   * present in the node.
   */
  void 
  insert(Vertex v,
         Filtration_value filtration_value)
  {
    typename Dictionary::iterator sh = members_.find(v);
    if(sh != members_.end() &&  sh->second.filtration() > filtration_value)
    { sh->second.assign_filtration(filtration_value);
      return; }
    if(sh == members_.end()) 
    {  members_.insert(std::pair< Vertex, Node >( v, Node(this,filtration_value) )); 
      return; }
  }

  typename Dictionary::iterator find(Vertex v)
  { return members_.find(v);  }
    
  /**********************/  
  
  Simplex_tree_siblings * oncles()
  {  return oncles_;  }
  
  Vertex parent()
  {  return parent_;  }
  
  Dictionary & members()
  { return members_; }
  
  size_t size() { return members_.size(); }


//private:
  Simplex_tree_siblings        * oncles_;
  Vertex                         parent_;
  Dictionary                     members_;
  
};

#endif // GUDHI_SIMPLEX_TREE_SIBLINGS
