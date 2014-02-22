/*
*  SimplexTree_iterators.h
*  Gudhi
*
*  Created by Clément Maria on 1/7/14.
*  Copyright 2014 INRIA. All rights reserved.
*
*/

#ifndef SIMPLEX_TREE_ITERATORS_H
#define SIMPLEX_TREE_ITERATORS_H

#include "boost/iterator/iterator_facade.hpp"

/** \brief Iterator over the vertices of a simplex
* in a SimplexTree.
*
* Forward iterator, 'value_type' is SimplexTree::Vertex.*/
template < class SimplexTree >
class Simplex_vertex_iterator 
: public boost::iterator_facade < Simplex_vertex_iterator < SimplexTree >
, typename SimplexTree::Vertex const
, boost::forward_traversal_tag
, typename SimplexTree::Vertex const
>
{
public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings       Siblings;
  typedef typename SimplexTree::Vertex         Vertex;

Simplex_vertex_iterator (SimplexTree * st) :   //any end() iterator
sib_(NULL), v_(st->null_vertex()) {}

Simplex_vertex_iterator( SimplexTree *  st,
  Simplex_handle sh) :
sib_(st->self_siblings(sh)),
v_(sh->first)    {}

private:
  friend class boost::iterator_core_access;

  bool equal (Simplex_vertex_iterator const &other) const
  { return sib_ == other.sib_ && v_ == other.v_; }

  Vertex const& dereference() const { return v_; }

  void increment () { v_ = sib_->parent(); sib_ = sib_->oncles();}

  Siblings *    sib_;
  Vertex        v_;
};

/*---------------------------------------------------------------------------*/
/** \brief Iterator over the simplices of the boundary of a
*  simplex.
*
* Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template < class SimplexTree >
class Boundary_simplex_iterator 
: public boost::iterator_facade < Boundary_simplex_iterator< SimplexTree >
, typename SimplexTree::Simplex_handle const
, boost::forward_traversal_tag
>
{
public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Vertex                Vertex;
  typedef typename SimplexTree::Siblings             Siblings;

// any end() iterator
  Boundary_simplex_iterator(SimplexTree * st) :
  last_(st->null_vertex()), sib_(NULL) {}

  Boundary_simplex_iterator(SimplexTree *  st,
    Simplex_handle sh) :
  suffix_(), st_(st)
  { 
    last_          = sh->first;
    Siblings * sib = st->self_siblings(sh);
    next_          = sib->parent();
sib_           = sib->oncles();       /** \todo check if NULL*/
    if(sib_ != NULL) { sh_ = sib_->find(next_); };
  }

private:
  friend class boost::iterator_core_access;
// valid when iterating along the SAME boundary.
  bool equal (Boundary_simplex_iterator const& other) const
  { return (sib_ == other.sib_ && last_ == other.last_);}

  Simplex_handle const& dereference () const  { return sh_; }

  void increment()
  { if(sib_ == NULL) { last_ = st_->null_vertex(); return; }

  Siblings * for_sib = sib_;
  for(typename std::vector< Vertex >::reverse_iterator rit = suffix_.rbegin();
    rit != suffix_.rend(); ++rit)
  {
    sh_ = for_sib->find(*rit);
    for_sib =  sh_->second.children();
  } 
sh_ = for_sib->find(last_); //sh_ points to the right simplex now
suffix_.push_back(next_);
next_ = sib_->parent();
sib_ = sib_->oncles();
}
Vertex                   last_   ; //last vertex of the simplex
Vertex                   next_   ; //next vertex to push in suffix_
std::vector< Vertex >    suffix_ ;
Siblings *               sib_    ; //where the next search will start from
Simplex_handle           sh_     ; //current Simplex_handle in the boundary
SimplexTree *            st_     ; //simplex containing the simplicial complex
};
/*---------------------------------------------------------------------------*/
/** \brief Iterator over the simplices of a simplicial complex.
*
* Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template < class SimplexTree >
class Complex_simplex_iterator 
: public boost::iterator_facade < Complex_simplex_iterator< SimplexTree >,
typename SimplexTree::Simplex_handle const, 
boost::forward_traversal_tag
>
{
public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings       Siblings;
  typedef typename SimplexTree::Vertex         Vertex;

//any end() iterator
  Complex_simplex_iterator() : st_(NULL) {}

  Complex_simplex_iterator(SimplexTree * st) :
  st_(st) 
  {
    if(st == NULL || st->root() == NULL || st->root()->members().empty())  { st_ = NULL; }
    else {
      sh_ = st->root()->members().begin();
      while(st->has_children(sh_)) 
        { sib_ = sh_->second.children();
          sh_ = sib_->members().begin();}
        }
      }
    private:
      friend class boost::iterator_core_access;

// valid when iterating along the SAME boundary.
      bool equal (Complex_simplex_iterator const& other) const
      { if(other.st_ == NULL) { return (st_ == NULL); }
      if(st_ == NULL) { return false; }
      return (&(sh_->second) == &(other.sh_->second));}

      Simplex_handle const& dereference () const { return sh_; }

// Depth first traversal.
      void increment ()
      { ++sh_;
        if(sh_ == sib_->members().end())
        {
if(sib_->oncles() == NULL) { st_ = NULL; return; } //reach the end
sh_ = sib_->oncles()->members().find(sib_->parent());
sib_ = sib_->oncles();    
return;  }
while(st_->has_children(sh_)) 
  { sib_ = sh_->second.children();
    sh_ = sib_->members().begin(); }
  }

  Simplex_handle   sh_;
  Siblings *               sib_;
  SimplexTree *            st_;
};

#endif // SIMPLEX_TREE_ITERATORS_H
