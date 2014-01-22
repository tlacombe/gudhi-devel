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

//#include <map>
#include "boost/container/flat_map.hpp"
#include "Filtered_simplex_tree_node.h"

/**
 * 
 */
class Simplex_tree_siblings {
public:
	
	typedef Filtered_simplex_tree_node														Node					;
	typedef int																										Vertex				;
	/** \brief Dictionary to store the nodes.
	 *
	 * Construct the relation Vertex -> Node, where the Vertex is the biggest 
	 * Vertex of the Simplex corresponding to the output Node.
	 * Must be ordered increasingly.
	 */
	typedef boost::container::flat_map<Vertex,Node>								Dictionary		;
//	typedef std::map<Vertex,Node>	Dictionary;
	
	typedef Dictionary::iterator																	Dictionary_it	;
	
	/**
	 * Default constructor
	 */
	Simplex_tree_siblings() :
	members_(Dictionary()),
	parent_(-1),
	oncles_(NULL)
	{};
	
	/**
	 * Construct with values
	 */
	Simplex_tree_siblings(Simplex_tree_siblings	* oncles,
												Vertex									parent ) :
	oncles_(oncles),
	parent_(parent),
	members_()
	{}
	
	
	/**
	 * 'members' is sorted and unique.
	 */
	Simplex_tree_siblings(Simplex_tree_siblings * oncles,
												Vertex									parent,
												std::vector< std::pair< Vertex, Node > > & members) :
	oncles_(oncles),
	parent_(parent),
	members_(boost::container::ordered_unique_range,members.begin(),members.end())
	{}
	
	/**
	 * Construct with initialized set of members
	 */
	Simplex_tree_siblings(Simplex_tree_siblings * oncles,
												Vertex									parent, 
												Dictionary						&	init_members) :
	oncles_(oncles),
	parent_(parent)
	{
		members_ = Dictionary(init_members);
		for(Dictionary::iterator map_it = members_.begin();
				map_it != members_.end(); map_it++)
		{	map_it->second.assign_children(this);	}
	}
		
	
	/**
	 * Destructor, calls recursively the destructor
	 * for all St_siblings
	 */
	~Simplex_tree_siblings()
	{
		for(Dictionary::iterator map_it = members_.begin();
				map_it != members_.end(); map_it++)
		{
			if(map_it->second.has_children(map_it->first)) 
			{
				delete map_it->second.children();
			}
		}
	};
	
	/**
	 * \brief Inserts a Node in the set of siblings nodes.
	 *
	 * If already present, assigns the minimal filtration value 
	 * between input filtration_value and the value already 
	 * present in the node.
	 */
	void 
	insert(Vertex v,
				 double filtration_value)
	{
		Dictionary_it sh = members_.find(v);
		if(sh != members_.end() &&  sh->second.filtration() > filtration_value)
		{ sh->second.assign_filtration(filtration_value);
			return; }
		if(sh == members_.end()) 
		{	members_.insert(std::pair< Vertex, Node >( v, Node(filtration_value,this) )); 
			return; }
	}
		
	
	Dictionary_it find(Vertex v)
	{ return members_.find(v);	}
	
	
	
	/**********************/	
	
	Simplex_tree_siblings * oncles()
	{	return oncles_;	}
	
	Vertex parent()
	{	return parent_;	}
	
	Dictionary & members()
	{ return members_; }
	
				
private:
	Simplex_tree_siblings				* oncles_;
	Vertex												parent_;
	Dictionary										members_;
	
};


#endif // GUDHI_SIMPLEX_TREE_SIBLINGS