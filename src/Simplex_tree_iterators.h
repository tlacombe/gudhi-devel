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

 /**
   * \brief Iterator over the vertices of a simplex
   * in a SimplexTree.
   *
   * 'value_type' must be Vertex.
   */
 template < class SimplexTree >
 class Simplex_vertex_iterator 
  : public boost::iterator_facade< Simplex_vertex_iterator< SimplexTree >,
                                   typename SimplexTree::Vertex,
                                   boost::forward_traversal_tag 
                                   >
  {
    public:
    typedef typename SimplexTree::Simplex_handle Simplex_handle;
    typedef typename SimplexTree::Siblings       Siblings;
    typedef typename SimplexTree::Vertex         Vertex;

    //any end() iterator
    Simplex_vertex_iterator() :
    sib_(NULL),
    v_(-1)
    {}

    Simplex_vertex_iterator(Simplex_handle sh) :
    sib_(sh->second.self_siblings(sh->first)),
    v_(sh->first)
    {}

    bool operator!= (const Simplex_vertex_iterator & other) const
    { if(v_ == other.v_ && sib_ == other.sib_) return false;
      return true; }

    Vertex operator* () const
    {return v_;}

    Simplex_vertex_iterator & operator++ ()
    { if(sib_ == NULL) 
        { v_ = -1;
          return *this; }
      v_ = sib_->parent();
      if(sib_->oncles() != NULL) { sib_ = sib_->oncles(); }
      else                       { sib_ = NULL; }    
      return *this;
    }

    private:
    Siblings *    sib_;
    Vertex        v_;
  };
/*------------------------------------*/
/**
   * \brief Range over the vertices of a simplex
   *
   * Methods .begin() and .end() return
   * a Simplex_vertex_iterator. 
   *
   * The order in which the Vertices are visited defines
   * the canonical orientation of the simplex.
   */
  template < class SimplexTree >
  class Simplex_vertex_range {
    public:
    typedef typename SimplexTree::Simplex_handle Simplex_handle;

    Simplex_vertex_range(Simplex_handle sh) :
    sh_(sh)
    {}

    Simplex_vertex_iterator < SimplexTree > begin ()
    { return Simplex_vertex_iterator < SimplexTree > (sh_); }

    Simplex_vertex_iterator < SimplexTree > end ()
    { return Simplex_vertex_iterator < SimplexTree > (); }
   
    private:
    Simplex_handle sh_;
  };

/*---------------------------------------------------------------------------*/
  /**
   * \brief Iterator over the simplices of the boundary of a
   * simplex.
   *
   * `value_type` must be `Simplex_handle`.
   */
  template < class SimplexTree >
  class Boundary_simplex_iterator 
    : public boost::iterator_facade < Boundary_simplex_iterator< SimplexTree >,
                                      typename SimplexTree::Simplex_handle,
                                      boost::forward_traversal_tag
                                    >
  {
    public:
    typedef typename SimplexTree::Simplex_handle Simplex_handle;
    typedef typename SimplexTree::Vertex         Vertex;
    typedef typename SimplexTree::Siblings       Siblings;

    Boundary_simplex_iterator() :
    sib_(NULL) {}
    
    Boundary_simplex_iterator(SimplexTree *  st,
                              Simplex_handle sh) :
    suffix_(),
    st_(st)
    {
        if(st == NULL) //end()
          { sib_ = NULL; }
        else 
    {
      last_          = sh->first;
      Siblings * sib = sh->second.self_siblings(last_);
      next_          = sib->parent();
      sib_           = sib->oncles();       /** \todo check if NULL*/
      sh_            = sib_->find(next_);
    }
  }
  /**
   * works only when compared with an end() iterator
   */
  bool operator!= (const Boundary_simplex_iterator & other) const
  { if(next_ < -1 ) return false; //indicates it's end()
    return true; } 

  Simplex_handle operator* () const
  { return sh_; }

  Boundary_simplex_iterator & operator++ ()
  { //stopping condition?
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
  Vertex                            last_   ; //last vertex of the simplex
  Vertex                            next_   ; //next vertex to push in suffix_
  std::vector< Vertex >             suffix_ ;
  Siblings *                        sib_    ; //where the next search will start from
  Simplex_handle                    sh_     ; //current Simplex_handle in the boundary
  SimplexTree *                     st_     ; //simplex containing the simplicial complex
};
/*-----------------------------------*/
 /**
   * \brief Range over the simplices in the boundary of a simplex.
   *
   * Methods .begin() and .end() return a Boundary_simplex_iterator.
   */
 template < class SimplexTree >
 class Boundary_simplex_range {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;

  Boundary_simplex_range(SimplexTree *    st,
                         Simplex_handle   sh) :
  st_(st),
  sh_(sh)
  {}

  Boundary_simplex_iterator < SimplexTree > begin()
  {return Boundary_simplex_iterator < SimplexTree > (st_,sh_);    }

  Boundary_simplex_iterator < SimplexTree > end()
  {return Boundary_simplex_iterator < SimplexTree > (NULL,sh_);}

 private:
  SimplexTree *       st_;
  Simplex_handle      sh_;
};

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
  /**
   * \brief Iterator over the simplices of a 
   * simplicial complex.
   *
   * 'value_type' must be Simplex_handle.
   */
template < class SimplexTree >
class Complex_simplex_iterator 
  : public boost::iterator_facade< Complex_simplex_iterator< SimplexTree >,
                                   typename SimplexTree::Simplex_handle, 
                                   boost::forward_traversal_tag
                                 >
{
  public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings       Siblings;
  typedef typename SimplexTree::Vertex         Vertex;
  typedef typename SimplexTree::Node           Node;

  Complex_simplex_iterator(Simplex_handle sh,
                          SimplexTree * st) :
  sh_(sh),
  st_(st)
  {}
  /**
   * Equal if points to the exact same place: the traversal is always
   * the same.
   */
  bool operator!= (const Complex_simplex_iterator & other) const
  { if(st_ == NULL) return false;
    //if( &(sh_->second) == &(other.sh->second) ) return false;
    return true;}

  Simplex_handle operator* () const
  { return sh_;}

  /** Depth first traversal.*/
  Complex_simplex_iterator & operator++ ()
  {
    if(sh_->second.has_children(sh_->first))
      { sh_ = sh_->second.children()->members().begin();
        return *this;}

    Siblings * sib = sh_->second.self_siblings(sh_->first);
    ++sh_;
    Vertex parent;

    while(sh_ == sib->members().end())
    {
      if(sib->oncles() == NULL)
      {
        parent = sib->parent();
        typename std::vector< Node >::iterator v_it = st_->root().begin()+parent+1;
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
  Simplex_handle        sh_;
  SimplexTree *         st_;
};
/*-----------------------------------*/
   /**
   * \brief Range over the simplices of a simplicial complex.
   *
   * Methods .begin() and .end() return
   * a Complex_simplex_iterator. 
   */
template < class SimplexTree >
class Complex_simplex_range {
  public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings       Siblings;
  typedef typename SimplexTree::Vertex         Vertex;
  typedef typename SimplexTree::Node           Node;

  Complex_simplex_range(SimplexTree * st) :
  st_(st) {}
  /**
   * \todo Problem here, should be able to have handles to vertices
   *  Must define a NULL handle too for a specific Simplex tree.
   */
  Complex_simplex_iterator < SimplexTree > begin()//traverse root until children
  {
    Vertex v = 0;
    typename std::vector< Node >::iterator v_it = st_->root().begin();
    for(;!(v_it->has_children(v)) && (v_it != st_->root().end()); ++v_it)
      {}

    if(v_it != st_->root().end()) 
    {
      return Complex_simplex_iterator < SimplexTree > 
                      (v_it->children()->members().begin(),
                       st_                                );
    }
    else {std::cerr << "Only Vertices in complex...";}
  }

  /** \todo implementation Complex_simplex_range.end() */
  Complex_simplex_iterator < SimplexTree > end();

 private:
  SimplexTree *    st_;
};

#endif // SIMPLEX_TREE_ITERATORS_H
