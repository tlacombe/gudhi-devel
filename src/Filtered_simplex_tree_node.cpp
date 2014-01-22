/*
 *  Filtered_simplex_tree_node.cpp
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#include "Simplex_tree_siblings.h"


Filtered_simplex_tree_node::
Filtered_simplex_tree_node() :
	children_(NULL),inter_node_(NULL),filtration_(0) {}

Filtered_simplex_tree_node::
Filtered_simplex_tree_node(double filtration,
							Simplex_tree_siblings *self_sib) :
children_(self_sib),
inter_node_(NULL),
filtration_(filtration) 
{}

Filtered_simplex_tree_node::
Filtered_simplex_tree_node(double filtration) :
children_(NULL),
inter_node_(NULL),
filtration_(filtration) 
{}
	
Filtered_simplex_tree_node::
Filtered_simplex_tree_node(const Filtered_simplex_tree_node &other) :
children_(other.children_),
inter_node_(other.inter_node_),
filtration_(other.filtration_)	
{}

Filtered_simplex_tree_node::
~Filtered_simplex_tree_node()	
{}


Simplex_tree_siblings *
Filtered_simplex_tree_node::self_siblings(Vertex label)
{
	Simplex_tree_siblings * next_sib = children_;
	if(next_sib == NULL) std::cerr << "Error in get_self_siblings \n";
	if(next_sib->parent() == label) return next_sib->oncles();
	else														return next_sib;
};

void 
Filtered_simplex_tree_node::list_of_vertices(std::vector<Vertex> &sigma,
											 Vertex							label)
{
	Simplex_tree_siblings *curr_sib = self_siblings(label);
	sigma.push_back(label);
	
	while(curr_sib != NULL)
	{
		sigma.push_back(curr_sib->parent());
		curr_sib = curr_sib->oncles();
	}
};

bool 
Filtered_simplex_tree_node::is_edge(Vertex label)
{
	if(self_siblings(label)->oncles() == NULL) return true;
	else return false;
};

bool 
Filtered_simplex_tree_node::has_children(Vertex label)
{
	if(children_ == NULL)							return false; //for root simplices
	if(children_->parent() == label)	return true;
	else															return false;
};

int 
Filtered_simplex_tree_node::dimension(unsigned int label)
{
	Simplex_tree_siblings *	curr_sib		= self_siblings(label);
	unsigned int						dim					= 0;
	
	while(curr_sib != NULL)
	{
		++dim;
		curr_sib = curr_sib->oncles();
	}
	return dim;
};