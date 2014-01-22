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
#include "Euclidean_rips_naive_geometry_traits.h"
#include "Filtered_simplex_tree_node.h"
#include "Simplex_tree_siblings.h"
#include "Flag_simplex_tree_iterators.h"



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
  typedef typename NeighborsGeometryTraits::FT              Filtration_value;
  typedef typename NeighborsGeometryTraits::Point           Point			;
  typedef typename NeighborsGeometryTraits::Point_range	    Point_range		;
  typedef typename NeighborsGeometryTraits::Point_iterator  Point_iterator	;
  typedef typename NeighborsGeometryTraits::Vertex          Vertex			; //-
  typedef typename NeighborsGeometryTraits::Vertex_iterator	Vertex_iterator	;
  typedef typename NeighborsGeometryTraits::Vertex_range    Vertex_range 	;
  typedef typename NeighborsGeometryTraits::Neighbor_vertex_range	
    												Neighbor_vertex_range	;
  typedef typename NeighborsGeometryTraits::Neighbor_vertex_iterator
    											Neighbor_vertex_iterator	;
  /** Node in the simplex tree.*/	
  typedef Filtered_simplex_tree_node						Node			;
  /** \brief Set of nodes sharing a same parent in the simplex tree. */
  typedef Simplex_tree_siblings								Siblings		;
  /** \brief Must be an ordered range. */
  typedef typename Siblings::Dictionary						Dictionary		;
  /** \todo Probably not correct, should be Dictionary::iterator. */
  typedef typename Siblings::Dictionary_it					Dictionary_it	;
	
  typedef typename Siblings::Dictionary_it					Simplex_handle	;
  /// @}


  /// \name Simplex Vertex traversal
  /// @{
  /**
  * \brief Range over the vertices of a simplex
  *
  * Methods .begin() and .end() return
  * a Simplex_vertex_iterator. 
  */
  class 											 		Simplex_vertex_range;
  /**
  * \brief Iterator over the vertices of a simplex
  *
  * 'value_type' must be Vertex.
  */
  class 													Simplex_vertex_iterator;
  /**
  * Returns a Simplex_vertex_range for the sequence of all
  * vertices of the simplex associated
  * to Simplex_handle sh.
  */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle sh)
  { return Simplex_vertex_range(sh); }

  void print(	Simplex_handle sh,	std::ostream& os = std::cout);
	
  /// @}






  /** \todo svr must be ordered */
  /*	Simplex_handle does_simplex_belong_to_complex(std::vector< Vertex > & svr)
	{	
	if(svr.begin() == svr.end()) std::cerr << "Empty simplex \n";
	 
	Simplex_vertex_iterator it = svr.begin();
	Simplex_handle sh = root_[*it];
	++it;
	for(;	it != svr.end(); ++it)
	{
	if(! sh->second.has_children(sh->first))
	return std::cerr << "Not here \n";       //must create a st.end() Simplex_handle...
	sh = sh->second.children()->find(*it);			// return some false if not here...
	}
	return sh;
	}
  */

  /** \todo svr must be ordered */



  Simplex_handle find(std::vector < Vertex > & s);
	 
	 
	 
  //	int simplex_dimension(Simplex & s) { return s.size()-1; }
	 
  /** \todo */
  /*	int complex_dimension()
	{ //...	
	}
  */



  /** \brief Default constructor.*/
  Flag_simplex_tree();

  /** Destructor.*/
  ~Flag_simplex_tree()		 { delete gt_; }


  /// \name Construction of the flag complex
  /// @{
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
	    int					dim_max,
	    Filtration_value	rho_max);
	
  /**
   * \brief Recursive expansion.
   */
  void siblings_expansion(Siblings * siblings, //must contain elements
			  int k);
	
  /**
   * Intersects Dictionary 1 [begin1;end1) with
   * Dictionary 2 [begin2,end2) and assign the
   * maximal possible value to the Nodes.
   */
  void intersection(std::vector< std::pair< Vertex, Node > > & intersection,
		    Dictionary_it		begin1,
		    Dictionary_it		end1,
		    Dictionary_it		begin2,
		    Dictionary_it		end2,
		    Filtration_value filtration);
  /// @}
	 
	 
  /// \name Acces methods
  /// @{
  /** Returns a pointer to the geometry traits.*/
  NeighborsGeometryTraits * gt ()			{ return gt_; }
  /** Returns the maximal threshold value.*/
  Filtration_value rho_max() 				{return rho_max_;}
  /** Returns the number of vertices in the complex.*/
  int 						nb_vertices() 		{return nb_vertices_; }
  /** \brief Returns the number of faces of the complex.*/
  int 						size_complex()	 { return size_cpx_; }
	 
  /** Returns a reference to the root nodes of the simplex tree.*/
  std::vector< Node > & root() {return root_; }
  /// @}
	 
	 
  /********** BOUNDARY SIMPLEX TRAVERSAL ***********/		
  /*************************************************/		
  /// \name Boundary Iterator
  /// @{
  /**
   * \brief Iterator on simplices belonging to a
   * boundary.
   *
   * `value_type` must be a `Simplex_handle`.
   */
  class Boundary_simplex_iterator;
	 
  /**
   * Range over the simplices in a boundary.
   * .begin() and .end() return Boundary_simplex_iterator
   * type.
   */
  class Boundary_simplex_range;
  /**
   *	\brief Range over all simplices of the boundary of a simplex, i.e.
   *   the set of codimension $1$ subsimplices of the Simplex.
   *
   * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
   * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
   * simplices of the boundary in the order: 
   * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
   *
   * We note that the alternate sum of the simplices given by the iterator
   * gives the chains corresponding to the boundary of the simplex.
   */
  Boundary_simplex_range boundary_simplex_range(Simplex_handle sh)
  { return Boundary_simplex_range(this,sh);}
  /// @}
  /*************************************************/		



  /*************************************************/		
  /// \name Filtration
  /// @{
  typedef std::vector< Simplex_handle >::iterator Filtration_simplex_iterator;

  class Filtration_simplex_range;


  std::vector< Simplex_handle > & filtration_vector()	{return filtration_vect_; }

  void initialize_filtration();


  /**
   * The use of stable_sort + the is_subface comparison
   * allows a nice filtration strategy order.
   */
  bool compare_simplices_fil(const Simplex_handle sh1,
			     const Simplex_handle sh2);

  /**
   * Returns true iff sh1 is a subface of sh2.
   */
  bool is_subface(Simplex_handle sh1, Simplex_handle sh2);
  /// @}





  /// \name Complex
  /// @{
  class Complex_simplex_iterator;
  class Complex_simplex_range;

  /**
   * Returns a range over all simplices in the complex
   * of dimension > 0
   */
  Complex_simplex_range complexe_simplex_range()
  {return Complex_simplex_range(this);}
  /// @}







 private:	
  /**
   * Maximum over 3 values.
   */
  Filtration_value maximum(Filtration_value a, 
			   Filtration_value b, 
			   Filtration_value c )
  {
    Filtration_value max = ( a < b ) ? b : a;
    return ( ( max < c ) ? c : max );
  }


  NeighborsGeometryTraits		*	gt_			;
  Filtration_value				rho_max_		;	
  int							nb_vertices_;
  int							size_cpx_		; //with or without vertices ?
  //	int							dimension_cpx_	; 
  std::vector< Node >				root_			;  //set of top nodes
	
  std::vector< Simplex_handle > 	filtration_vect_;

  Simplex_handle					NULL_sh_		;
};






template<class N>
class Flag_simplex_tree<N>::Simplex_vertex_iterator {
 public:
	
  //any end() iterator
 Simplex_vertex_iterator() :
  sib_(NULL),
    v_(-1)
      {}
	
 Simplex_vertex_iterator(Simplex_handle sh) :
  sib_(sh->second.self_siblings(sh->first)),
    v_(sh->first)
      {}
	
  /**
   * 
   */
  bool operator!= (const Simplex_vertex_iterator & other) const
  {
    if(v_ == other.v_ && sib_ == other.sib_) return false;
    return true;
  }

  Vertex operator* () const
  {return v_;}

  const Simplex_vertex_iterator & operator++ ()
  {
    if(sib_ == NULL) 
      {
	v_ = -1;
	return *this;
      }
    v_ = sib_->parent();
    if(sib_->oncles() != NULL)	{	sib_ = sib_->oncles();}
    else						{	sib_ = NULL	;}	
    return *this;
  }

 private:
  Siblings *	sib_;
  Vertex			v_;
};

template<class N>
class Flag_simplex_tree<N>::Simplex_vertex_range {
 public:

  typedef typename Simplex_tree_siblings::Dictionary_it Simplex_handle;

 Simplex_vertex_range(Simplex_handle sh) :
  sh_(sh)
  {}

  Simplex_vertex_iterator begin ()
  {return Simplex_vertex_iterator(sh_);}

  Simplex_vertex_iterator end ()
  { return Simplex_vertex_iterator(); }
 private:
  Simplex_handle sh_;
};






template<class N>
class Flag_simplex_tree<N>::Boundary_simplex_iterator {
 public:
 Boundary_simplex_iterator() :
  sib_(NULL)
    {}

  /**
   */
 Boundary_simplex_iterator(	Flag_simplex_tree * st,
				Simplex_handle sh) :
  suffix_(),
    st_(st)
    {
      if(st == NULL) //end()
	{	sib_ = NULL;	}
      else 
	{
	  last_			= sh->first;
	  Siblings * sib	= sh->second.self_siblings(last_);
	  next_			= sib->parent();
	  sib_			= sib->oncles();       /** \todo check if NULL*/
	  sh_				= sib_->find(next_);
	}
    }
  /**
   * works only when compared with an end() iterator
   */
  bool operator!= (const Boundary_simplex_iterator & other) const
  {
    if(next_ < -1 ) return false; //indicates it's end()
    return true;
  }	

  Simplex_handle operator* () const
  { return sh_; }

  const Boundary_simplex_iterator & operator++ ()
  {
    //stopping condition?
    if(sib_ == NULL) {next_ = -2; return *this;}
	
    // search the new Simplex_handle
    Siblings * for_sib = sib_;
    for(typename std::vector< Vertex >::reverse_iterator rit = suffix_.rbegin();
	rit != suffix_.rend(); ++rit)
      {
	sh_ = for_sib->find(*rit);
	for_sib =  sh_->second.children();
      } 
    sh_ = for_sib->find(last_); //sh_ points to the right simplex now
	
    if(sib_->oncles() != NULL) //middle of the tree
      {	
	suffix_.push_back(next_);
	next_ = sib_->parent();
	sib_ = sib_->oncles();
      }
    else // special case sib_->oncles() == root
      { 
	if(next_ == -1) {sib_ = NULL; return *this;}
	sib_ = st_->root()[next_].children(); //trick
	next_ = -1;          // ++ then ++ and it's .end()
      }
    return *this;
  }

 private:
  Vertex							last_	; //last vertex of the simplex
  Vertex							next_   ; //next vertex to push in suffix_
  std::vector< Vertex >			suffix_	;
  Siblings					* 	sib_	; //where the next search will start from
  Simplex_handle 					sh_		; //current Simplex_handle in the boundary
  Flag_simplex_tree 			* 	st_		; //simplex containing the simplicial complex
};

/**
 * Range over the simplices in a boundary.
 * .begin() and .end() return Boundary_simplex_iterator
 * type.
 */
template<class N>
class Flag_simplex_tree<N>::Boundary_simplex_range {
 public:
 Boundary_simplex_range(Flag_simplex_tree *	st,
			Simplex_handle				sh) :
  st_(st),
    sh_(sh)
    {}

  Boundary_simplex_iterator begin()
  {return Boundary_simplex_iterator(st_,sh_);	}

  Boundary_simplex_iterator end()
  {return Boundary_simplex_iterator(NULL,sh_);}

 private:
  Flag_simplex_tree		* 	st_;
  Simplex_handle				sh_;
};







template<class N>
class Flag_simplex_tree<N>::Filtration_simplex_range {
 public:
 Filtration_simplex_range(Flag_simplex_tree * st) :
  st_(st)
  {}

  /**
   * Initializes the filtration_vect too
   */
  Filtration_simplex_iterator begin()
  {
    if(st_->filtration_vector().size() > 0)
      { //initialize the filtration vector of st
	st_->initialize_filtration();
      }
    return st_->filtration_vector().begin();
  }

  Filtration_simplex_iterator end()
  {st_->filtration_vector().end();}

 private:
  Flag_simplex_tree  	*	st_;
};









template<class N>
class Flag_simplex_tree<N>::Complex_simplex_iterator {
 public:
 Complex_simplex_iterator(Simplex_handle sh,
			  Flag_simplex_tree * st) :
  sh_(sh),
    st_(st)
    {}

  /**
   * Equal if point to the exact same place: the traversal is always
   * the same.
   */
  bool operator!= (const Complex_simplex_iterator & other) const
  {
    if(st_ == NULL) return false;
    if( &(sh_->second) == &(other.sh->second) ) return false;
    return true;
  }

  Simplex_handle operator* () const
  { return sh_;}

  /**
   * Depth first traversal.
   */
  const Complex_simplex_iterator & operator++ ()
  {
    if(sh_->second.has_children(sh_->first))
      {
	sh_ = sh_->second.children()->members().begin();
	return *this;
      }

    Siblings * sib = sh_->second.self_siblings(sh_->first);
    ++sh_;

    typename Flag_simplex_tree<N>::Vertex parent;

    while(sh_ == sib->members().end())
      {
	if(sib->oncles() == NULL)
	  {
	    parent = sib->parent();
	    std::vector< Node >::iterator v_it = st_->root().begin()+parent+1;
	    ++parent;
	    while(v_it != st_->root().end() && !(v_it->has_children(parent)))
	      {++v_it; ++parent;}
	    if(v_it == st_->root().end()) {st_ = NULL;}
	    else {sh_ = v_it->children()->members().begin();}
	    return *this;
	  }
		
	else
	  {
	    parent = sib->parent();
	    sib = sib->oncles();
	    sh_ = sib->find(parent);
	    ++sh_;
	  }
      }
    return *this;
  }

 private:
  Simplex_handle		sh_;
  Flag_simplex_tree  * st_;
};

template<class N>
class Flag_simplex_tree<N>::Complex_simplex_range {
 public:
 Complex_simplex_range(Flag_simplex_tree * st) :
  st_(st)
  {}

  /**
   * \todo Problem here, should be able to have handles to vertices
   * Must define a NULL handle too for a specific Simplex tree.
   */
  Complex_simplex_iterator begin()//traverse root until children
  {
    Vertex v = 0;
    std::vector< Node >::iterator v_it = st_->root().begin();
    for(;!(v_it->has_children(v)) && (v_it != st_->root().end()); ++v_it)
      {}

    if(v_it != st_->root().end()) {return Complex_simplex_iterator(v_it->children()->members().begin(),
								   st_);}
    else {std::cerr << "Only Vertices in complex...";}
  }

  Complex_simplex_iterator end();

 private:
  Flag_simplex_tree * st_;
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


#include "Flag_simplex_tree.hpp"


#endif // GUDHI_FLAG_SIMPLEX_TREE_H
