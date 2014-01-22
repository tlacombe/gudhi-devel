/*
 *  Euclidean_rips_naive_geometry_traits.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/8/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_EUCLIDEAN_NAIVE_GEOMETRY_TRAITS_H
#define GUDHI_EUCLIDEAN_NAIVE_GEOMETRY_TRAITS_H

#include <set>
#include <cmath>

/**
 * \brief Represents Points in a Euclidean space with
 * nearest neighbors queries within a radius of a Point.
 *
 * The computation of nearest neighbors is naive, and 
 * traverses all Points in the metric space to find the one
 * within a given radius.
 * In particular, it does not store any additional data 
 * structure apart from the set of points.
 *
 * \todo Make it lazy.
 *
 * \implements NeighborsGeometryTraits
 */
class Euclidean_rips_naive_geometry_traits {
	public:

	/** \brief Distance type.*/
	typedef double														FT 	;
	typedef std::vector<double>								Point			;
	typedef std::vector< Point >							Point_range		;
	typedef std::vector< Point >::iterator					Point_iterator	;

	/** \brief Vertex type.*/
	typedef int												Vertex			;

	// Vertex_iterator type -----------------------------------------------------
	class Vertex_iterator {
	public:
		Vertex_iterator(Vertex idx) :
		v_(idx)
		{}
	
		bool operator!= (const Vertex_iterator& other) const
		{	return v_ != other.v_;}
		
		Vertex operator* ()
		{return v_;}
		
		Vertex_iterator & operator++ ()
		{	++v_;
			return *this;}
		
	private:
		Vertex	v_;
	}; //------------------------------------------------------------------------
	
	// --------------------------------------------------------------------------
	class Vertex_range {
	public:
		Vertex_range(int nb_vertices) :
		max_v_(nb_vertices)
		{}
		
		Vertex_iterator begin ()
		{ return Vertex_iterator(0); }
		
		Vertex_iterator	end ()
		{ return Vertex_iterator(max_v_); }
		
	private:
		Vertex			max_v_;
	}; //------------------------------------------------------------------------
		
	typedef std::vector< Vertex >::iterator						Neighbor_vertex_iterator;
	
	class Neighbor_vertex_range {
	public:
		Neighbor_vertex_range(Euclidean_rips_naive_geometry_traits * gt,
													 Vertex v,
													 FT dist_max) :
		v_(v),
		dist_max_(dist_max),
		gt_(gt),
		neighbors_(std::vector< Vertex >())
		{}
		
		Neighbor_vertex_iterator begin ()
		{
			double dist_squared = dist_max_ * dist_max_;
			Point p = gt_->vertex_to_point(v_);
			Vertex_iterator v_end = gt_->vertex_range().end();
			for(Vertex_iterator v_it = gt_->vertex_range().begin();
					v_it != v_end; ++v_it)
			{	//predicate
				if(v_ != *v_it && 
					 gt_->squared_distance(p,gt_->vertex_to_point(*v_it)) <= dist_squared)
				{ neighbors_.push_back (*v_it); }
			}
			return neighbors_.begin();
		}
		
		Neighbor_vertex_iterator end ()
		{ return neighbors_.end(); }
	
	private:
		Vertex																	v_;
		FT																			dist_max_;
		std::vector< Vertex >										neighbors_;
		Euclidean_rips_naive_geometry_traits	* gt_;
	}; //------------------------------------------------------------------------
	
	
	
	/**
	 * Default constructor
	 */
	Euclidean_rips_naive_geometry_traits() :
	dimension_(-1),
	nb_elements_(0),
	point_range_ptr_(NULL),
	vertex_range_(Vertex_range(0))
//	dist_max_(0)
	{}
	
	/**
	 * Initialize the trait.
	 *
	 * For example, construct a kd-tree. Here, copies a
	 * pointer to the set of Points.
	 */
	void 
	init(Point_range &points)
	{	
		nb_elements_ = points.size();
		point_range_ptr_ = &points;
		vertex_range_ = Vertex_range((int)points.size());
		if(! points.empty()) dimension_ = points.begin()->size();
		else dimension_ = -1;
	}
		
	
	
	/**
	 * Returns the point corresponding to a given
	 * Vertex. This allows the correspondance
	 * between combinatorial vertices and geometric
	 * points.
	 */
	Point &vertex_to_point(Vertex v)
	{	return (*point_range_ptr_)[v]; }
	
	/**
	 * \todo Necessary?
	 */
	Vertex point_to_vertex(Point & p);
	

	
	
	
	/**
	 * \brief Euclidean distance squared.
	 *
	 * Faster to compute than Euclidean distance and
	 * enough for comparison.
	 */
	FT squared_distance (Point &p1,
											 Point &p2)
	{	double squared_distance = 0.;
		Point::iterator p1_it = p1.begin();
		Point::iterator p2_it = p2.begin();
		for(; p1_it != p1.end(); p1_it++, p2_it++)
		{	squared_distance += (*p1_it - *p2_it) * (*p1_it - *p2_it); }
		return squared_distance; }
	
	/** \brief Euclidean distance.*/
	FT distance ( Point &p1, Point &p2 )
	{	return std::sqrt(squared_distance(p1,p2)); }
	
	/** \brief Distance between the Points corresponding to two Vertices.*/
	FT distance (Vertex u, Vertex v)
	{ return distance(vertex_to_point(u),vertex_to_point(v)); }
		
	
	/** \brief Returns the dimension of the Euclidean space.*/
	int dimension()	{ return dimension_; }
	
	/** \brief Returns a range over all Points.
	 *
	 * \todo not right, returns a ref.
	 */
	Point_range & point_range()	{ return * point_range_ptr_; }
	
	/** \brief Returns a range of all Vertices.*/
	Vertex_range vertex_range()	{ return vertex_range_; }

	
	int nb_elements()
	{ return nb_elements_; }
				
	private :
	/** Dimension of the embedding Euclidean space.*/
	int							dimension_;    
	int							nb_elements_;
	/** Range of points contained in the space. */
	Point_range		*	point_range_ptr_;
	Vertex_range    vertex_range_;
//	FT							dist_max_;
	
};

#endif // GUDHI_EUCLIDEAN_NAIVE_GEOMETRY_TRAITS_H