/*
 *  Simplex_tree_node.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_FILTERED_SIMPLEX_TREE_NODE_H
#define GUDHI_FILTERED_SIMPLEX_TREE_NODE_H

#include <vector>
#include <iostream>

class Simplex_tree_siblings;

/**
 * \brief Node of a simplex tree with filtration value
 * and simplex data.
 *
 * It contains its own filtration value.
 */
class Filtered_simplex_tree_node {
	public : 
	typedef int								Vertex;
	/**
	 * Default constructor.
	 */
	Filtered_simplex_tree_node() :
	children_(NULL),
	inter_node_(NULL),
	filtration_(0)
	{}
	
	/**
	 * Constructor with values.
	 */
	Filtered_simplex_tree_node(double filtration,
														 Simplex_tree_siblings *self_sib) :
	children_(self_sib),
	inter_node_(NULL),
	filtration_(filtration)
	{}
	
	Filtered_simplex_tree_node(double filtration) :
	children_(NULL),
	inter_node_(NULL),
	filtration_(filtration)
	{}
		
	/** Necessary
	 * Copy constructor
	 */
	Filtered_simplex_tree_node(const Filtered_simplex_tree_node &other) :
	children_(other.children_),
	inter_node_(other.inter_node_),
	filtration_(other.filtration_)
	{}
	
	/**
	 * Destructor
	 */
	~Filtered_simplex_tree_node()
	{}
	
	/**
	 * When in a simplex tree, returns a pointer
	 * to the Simplex_tree_siblings containing the node.
	 *
	 * Vertex label is the biggest Vertex in the simplex
	 * represented by the node.
	 */
	Simplex_tree_siblings *get_self_siblings(Vertex label);
	
	/**
	 * Return true if the node has children,
	 * false otherwise.
	 */
	bool has_children(Vertex label);
		
	/**
	 * Assign a children to the node
	 */
	void assign_children(Simplex_tree_siblings *children)
	{	children_ = children;	}
	
	/**
	 *
	 */
	void assign_filtration(double filtration_value)
	{	filtration_ = filtration_value;	}
	
	/**
	 * Return the dimension of the simplex corresponding
	 * to the node.
	 */
	int dimension(unsigned int label);
	
	/**
	 * Fill sigma with the list of vertices
	 * of the simplex corresponding to the node
	 */		
	void list_of_vertices(std::vector<Vertex> &sigma,
												Vertex							label);

	
	/**
	 * Return true is the simplex corresponding 
	 * to the node is an edge, false otherwise.
	 */
	bool is_edge(Vertex label);
	
	double filtration()
	{ return filtration_; }

	Simplex_tree_siblings * children()
	{ return children_; }
	
	void * inter_node()
	{ return inter_node_; }
	
	/***************************/
	private :	
	Simplex_tree_siblings		* children_;
	//S_inter_Node					* inter_node_;
	void										* inter_node_;
	double										filtration_; //value in the filtration
	
};


#endif // GUDHI_FILTERED_SIMPLEX_TREE_NODE_H

