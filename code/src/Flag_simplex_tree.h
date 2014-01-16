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
#include "Filtered_simplex_tree_node.h"

#include "Euclidean_rips_naive_geometry_traits.h"


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
	typedef typename NeighborsGeometryTraits::FT							Filtration_value	;
	typedef typename NeighborsGeometryTraits::Point						Point							;
	typedef typename NeighborsGeometryTraits::Point_range			Point_range				;
	typedef typename NeighborsGeometryTraits::Point_iterator	Point_iterator		;
	typedef typename NeighborsGeometryTraits::Vertex					Vertex						; //-
	typedef typename NeighborsGeometryTraits::Vertex_iterator	Vertex_iterator		;
	typedef typename NeighborsGeometryTraits::Vertex_range		Vertex_range			;
	typedef typename NeighborsGeometryTraits::Neighbor_vertex_range	
																										Neighbor_vertex_range			;
	typedef typename NeighborsGeometryTraits::Neighbor_vertex_iterator
																										Neighbor_vertex_iterator	;
	/** Node in the simplex tree.*/	
	typedef Filtered_simplex_tree_node												Node							;
	/** \brief Set of nodes sharing a same parent in the simplex tree. */
	typedef Simplex_tree_siblings															Siblings					;
	/** \brief Must be an ordered range. */
	typedef typename Siblings::Dictionary											Dictionary				;
	/** \todo Probably not correct, should be Dictionary::iterator. */
	typedef typename Siblings::Dictionary_it									Dictionary_it			;
	
	
	
	typedef typename	Siblings::Dictionary_it						Simplex_handle					;
//	typedef std::vector< Vertex >											Simplex									;
	typedef						std::vector< Vertex >							Simplex_vertex_range		;
	typedef typename	std::vector< Vertex >::iterator		Simplex_vertex_iterator	;
	
	/// @}

	
	/** \todo */
/*	bool does_simplex_belong_to_complex(Simplex & s)
	{	//...
	}
	*/
//	int simplex_dimension(Simplex & s) { return s.size()-1; }
	
	/** \todo */
/*	int complex_dimension()
	{ //...	
	}
	*/
	
	
	
	/** \brief Default constructor.*/
	Flag_simplex_tree() :
	gt_(new NeighborsGeometryTraits()),
	rho_max_(0),
	nb_vertices_(0),
	size_cpx_(0),
	root_()
	{};
	
	/** Destructor.*/
	~Flag_simplex_tree()
	{ delete gt_; }
		
	/**
	 * \brief Construct the flag complex.
	 *
	 * First introduces all edges given by Neighbor_vertex_range
	 * defined in NeighborsGeometryTraits: it produces, for a given
	 * vertex, a range allowing to iterate over its neighbors in 
	 * the 1-skeleton.
	 * Then, realizes the expansion of the graph until dimension
	 * dim_max.
	 */
	void init(//Point_range_sc	&	point_range,
						int								dim_max,
						Filtration_value	rho_max)
	{
		rho_max_ = rho_max;
		nb_vertices_ = gt_->nb_elements();
		// Insert all edges
		root_ = std::vector< Node >( gt_->nb_elements(), Node() );	
		for(Vertex_iterator v_it = gt_->vertex_range().begin();
				v_it != gt_->vertex_range().end(); ++v_it)
		{
			Neighbor_vertex_range n_range(gt_,*v_it,rho_max);
			for(Neighbor_vertex_iterator n_it = n_range.begin();
					n_it != n_range.end(); ++n_it)
			{
				if(*v_it < *n_it) 
				{
					if(! root_[*v_it].has_children(*v_it)) 
					{	root_[*v_it].assign_children(new Siblings(NULL,*v_it)); }
					root_[*v_it].children()->insert(*n_it,gt_->distance(*v_it,*n_it));
				}
			}
		}
		// Update size of the complex
		size_cpx_ += root_.size();
		int v = 0;
		for(std::vector< Node >::iterator r_it = root_.begin();
				r_it != root_.end(); ++r_it,++v)
		{ if(r_it->has_children(v)) {size_cpx_ += r_it->children()->members().size();}	}
	
		// Expansion
		clock_t start = clock();
		int curr_vertex = 0;
		for(std::vector< Node >::iterator root_it = root_.begin();
				root_it != root_.end(); ++root_it, ++curr_vertex)
		{
			if(root_it->has_children(curr_vertex)) 
			{ siblings_expansion(root_it->children(), dim_max-1); }
		}
		
		clock_t end = clock();
		std::cout << "Computational time for Rips construction = " << 
		(double)(end - start)/(double)CLOCKS_PER_SEC << std::endl;
	}
	
	/**
	 * \brief Recursive expansion.
	 */
	void siblings_expansion(Siblings * siblings, //must contain elements
													int k)
	{
		//	if (k==0 || members_.empty()) return;
		if(k == 0) return;
		
		Dictionary_it next = siblings->members().begin(); ++next;
		
		static std::vector< std::pair<Vertex , Node> > inter;

		for(Dictionary_it s_h = siblings->members().begin();
				s_h != siblings->members().end(); ++s_h,++next)
		{
			if(root_[s_h->first].has_children(s_h->first))
			{
				int size_intersection = intersection(inter,  //output intersection
																						 next,					//begin
																						 siblings->members().end(),//end
																						 root_[s_h->first].children()->members().begin(),
																						 root_[s_h->first].children()->members().end(),
																						 s_h->second.filtration());
				if(size_intersection != 0)
				{
					size_cpx_ += size_intersection;
					Siblings * new_sib = new Siblings(siblings,					 //oncles
																						s_h->first,		 //parent
																						size_intersection,
																						inter);					//boost::container::ordered_unique_range_t
					inter.clear();
					s_h->second.assign_children(new_sib);
					siblings_expansion(new_sib,k-1);
				}
				else {	inter.clear();	}
			}
		}
	}
	
		
	/**
	 * Intersects Dictionary 1 [begin1;end1) with
	 * Dictionary 2 [begin2,end2) and assign the
	 * maximal possible value to the Nodes.
	 */
	int intersection(std::vector< std::pair< Vertex, Node > > & intersection,
									 Dictionary_it		begin1,
									 Dictionary_it		end1,
									 Dictionary_it		begin2,
									 Dictionary_it		end2,
									 Filtration_value filtration)
	{
		if(begin1 == end1 || begin2 == end2) return 0;
		int size_intersection = 0;
		while( true )
		{
			if( begin1->first == begin2->first )
			{
				++size_intersection;
				intersection.push_back(std::pair< Vertex, Node >(begin1->first,
																												 Node(maximum(begin1->second.filtration(),
																																			begin2->second.filtration(),
																																			filtration))));
				++begin1;
				++begin2;
				if( begin1 == end1 || begin2 == end2 ) return size_intersection;
			}
			else { 
				if( begin1->first < begin2->first ) 
				{
					++begin1;
					if(begin1 == end1) return size_intersection;
				}
				else {
					++begin2;
					if(begin2 == end2) return size_intersection;
				}
			}
		}
	}
			
	/// \name Acces methods
	/// @{
	/** Returns a pointer to the geometry traits.*/
	NeighborsGeometryTraits * gt ()	{ return gt_; }
	/** Returns the maximal threshold value.*/
	Filtration_value rho_max() {return rho_max_;}
	/** Returns the number of vertices in the complex.*/
	int nb_vertices() {return nb_vertices_; }
	/** \brief Returns the number of faces of the complex.*/
	int size_complex() { return size_cpx_; }
	
	/** Returns a reference to the root nodes of the simplex tree.*/
	std::vector< Node > & root() {return root_; }
	/// @}
	
	
	
private:	
	/**
	 * Maximum over 3 values.
	 */
	Filtration_value maximum( Filtration_value a, 
													 Filtration_value b, 
													 Filtration_value c )
	{
		Filtration_value max = ( a < b ) ? b : a;
		return ( ( max < c ) ? c : max );
	}
	
	
	private :
	NeighborsGeometryTraits		*	gt_						;
	Filtration_value						rho_max_			;	
	int													nb_vertices_	;
	int													size_cpx_			; //
//	int													dimension_cpx_; 
	std::vector< Node >					root_					;	       //set of top nodes
	
};









std::ostream& operator<<(std::ostream& os, Simplex_tree_siblings & obj)
{
	os << "--Oncles: @ " << (long int)(obj.oncles()) << "\n";
	os << "--Parent:   " << obj.parent() << "\n";
	os << "Siblings: @ " << (long int)(&obj) << "\n";
	for(Simplex_tree_siblings::Dictionary_it sh = obj.members().begin();
			sh != obj.members().end(); ++sh)
	{	os << "[" << sh->first << ":" << sh->second.filtration() <<"] ";	}
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
//	os << "Dimension  = " << obj.dimension_cpx_ << std::endl;
	os << "rho_max = " << obj.rho_max() << std::endl;
	os << "nb_V    = " << obj.nb_vertices() << std::endl;
	os << std::endl;
	
	int v = 0;
	os << "@ 0000000000:   ";
	for(std::vector< Filtered_simplex_tree_node >::iterator it = obj.root().begin();
			it != obj.root().end(); ++it, ++v)
	{	os << v << " ";}
	os << std::endl;
	
	v = 0;
	for(std::vector< Filtered_simplex_tree_node >::iterator it = obj.root().begin();
			it != obj.root().end(); ++it,++v)
	{
		if(it->has_children(v)) os << *(it->children());	
	}
	return os;
};	





#endif // GUDHI_FLAG_SIMPLEX_TREE_H