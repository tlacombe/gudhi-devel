/*
 *  Hasse_diagram.h
 *  Gudhi
 *
 *  Created by Clément Maria on 4/2/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_HASSE_DIAGRAM_H
#define GUDHI_HASSE_DIAGRAM_H

#include "boost/iterator/counting_iterator.hpp"


template < typename FiltrationValueType = double
         , typename SimplexKeyType      = int    //must be a signed integer type
         , typename VertexHandleType    = int    //must be a signed interger type, int convertible to it
         >
class Hasse_complex 
{
public:

template < class HasseCpx >
struct Hasse_simplex {
  //Complex_ds must verify that cpx->key(sh) is the order of sh in the filtration
    template< class Complex_ds >
    Hasse_simplex ( Complex_ds *                        cpx
                  , typename Complex_ds::Simplex_handle sh )
    : key_(cpx->key(sh))
    , filtration_(cpx->filtration(sh))
    , boundary_()
    {
      boundary_.reserve(cpx->dimension(sh)+1);
      for( auto b_sh : cpx->boundary_simplex_range(sh) )
      { boundary_.push_back( cpx->key(b_sh) ); }
    }

    typename HasseCpx::Simplex_key                 key_;
    typename HasseCpx::Filtration_value            filtration_;
    std::vector<typename HasseCpx::Simplex_handle> boundary_;
  };

  typedef Hasse_simplex<Hasse_complex>         Hasse_simp;
  typedef FiltrationValueType                  Filtration_value;
  typedef SimplexKeyType                       Simplex_key;
  typedef int                                  Simplex_handle; //index in vector complex_

  typedef boost::counting_iterator< Simplex_handle >          Filtration_simplex_iterator;
  typedef boost::iterator_range<Filtration_simplex_iterator>  Filtration_simplex_range;  

  typedef typename std::vector< Simplex_handle >::iterator    Boundary_simplex_iterator;
  typedef boost::iterator_range<Boundary_simplex_iterator>    Boundary_simplex_range;  

  typedef typename std::vector< Simplex_handle >::iterator    Skeleton_simplex_iterator;
  typedef boost::iterator_range< Skeleton_simplex_iterator >  Skeleton_simplex_range;


/** \todo NOT VALID Hasse_complex skeleton_simplex_range(...) */
  Skeleton_simplex_range skeleton_simplex_range( int dim = 0 ) {
    if(dim != 0) { std::cerr << "Dimension must be 0 \n"; }
    return Skeleton_simplex_range(vertices_.end(),vertices_.end());
  }

  template < class Complex_ds >
  Hasse_complex(Complex_ds * cpx)
  : complex_()
  , vertices_()
  , threshold_(cpx->filtration())
  , num_vertices_(0)
  , dim_max_(cpx->dimension())
  {
    complex_.reserve(cpx->num_simplices());
    int idx = 0;
    for(auto cpx_sh : cpx->filtration_simplex_range())
    { 
      complex_.push_back(Hasse_simp(cpx,cpx_sh)); 
      if(dimension(idx) == 0) { vertices_.push_back(idx); } 
      ++idx; 
    }
  }

  size_t num_simplices() { return complex_.size(); }

  Filtration_simplex_range filtration_simplex_range() 
  { return Filtration_simplex_range( Filtration_simplex_iterator(0)
                                   , Filtration_simplex_iterator(complex_.size()) ); }

  Simplex_key key( Simplex_handle sh ) { return complex_[sh].key_; }

  Simplex_key null_key() { return -1; }

  Simplex_handle simplex( Simplex_key key ) 
  {
    if(key == null_key()) return null_simplex();
    return key;
  }

  Simplex_handle null_simplex() { return -1; }

  Filtration_value filtration( Simplex_handle sh ) {
    if( sh == null_simplex() ) { return filtration(); }
    return complex_[sh].filtration_;
  }

  Filtration_value filtration() { return threshold_; }

  int dimension ( Simplex_handle sh ) { return complex_[sh].boundary_.size()-1; }
  int dimension () { return dim_max_; }

  std::pair<Simplex_handle,Simplex_handle> endpoints( Simplex_handle sh ) 
  { return std::pair<Simplex_handle,Simplex_handle>( complex_[sh].boundary_[0]
                                                   , complex_[sh].boundary_[1] ) ;}

  void assign_key( Simplex_handle sh, Simplex_key key) { complex_[sh].key_ = key; }

  Boundary_simplex_range boundary_simplex_range ( Simplex_handle sh ) 
  { return Boundary_simplex_range( complex_[sh].boundary_.begin()
                                 , complex_[sh].boundary_.end() ); }


  std::vector< Hasse_simp >   complex_;
  std::vector<Simplex_handle> vertices_;
  Filtration_value            threshold_;
  size_t                      num_vertices_;
  int                         dim_max_;
};

#endif // GUDHI_HASSE_DIAGRAM_H

//concept PersistentCohomologySimplicialComplex
// template < typename FiltrationValueType = double
//          , typename SimplexKeyType      = int    //must be a signed integer type
//          , typename VertexHandleType    = int    //must be a signed interger type, int convertible to it
// //         , bool ContiguousVertexHandles = true   //true is Vertex_handles are exactly the set [0;n)
//          >
// class Hasse_diagram {

// typedef    ...    Filtration_simplex_range;  
// typedef    ...    Simplex_handle;
// typedef    ...    Skeleton_simplex_range;
// typedef    ...    Boundary_simplex_range;

// typedef FiltrationValueType Filtration_value;
// typedef SimplexKeyType Simplex_key;

// size_t num_simplices() {}

// Filtration_simplex_range filtration_simplex_range() {}

// Simplex_key key( Simplex_handle sh ) {}

// Simplex_key null_key() {}

// Simplex_handle simplex( Simplex_key key ) {}

// Simplex_handle null_simplex() {}

// Filtration_value filtration( Simplex_handle sh ) {
//   if( sh == null_simplex() ) { return filtration(); }
// }
// Filtration_value filtration() {}

// int dimension ( Simplex_handle ) {}

// Skeleton_simplex_range skeleton_simplex_range( int dim = 0 ) {}

// std::pair<Simplex_handle,Simplex_handle> endpoints( Simplex_handle sh ) {}

// void assign_key( Simplex_handle sh, Simplex_key key) {}

// Boundary_simplex_range boundary_simplex_range ( Simplex_handle sh ) {}

// };
