/*
 *  Persistence_cohomology.h
 *  Gudhi
 *
 *  Created by Clément Maria on 10/19/12.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

#ifndef _PERSISTENCECOMPUTATION_SIMPLEXTREE_
#define _PERSISTENCECOMPUTATION_SIMPLEXTREE_

#include <boost/progress.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/intrusive/list.hpp>

#include "Column_list.h"
#include "Multi_field.h"

#include <boost/pool/object_pool.hpp>

/** \brief Computes the persistent cohomology of a filtered simplicial complex.
*
* The computation is implemented with a Compressed Annotation Matrix.
* FilteredSimplicialComplexDS::Simplex_key must be int.
*
* \todo Filter for simplices (dimension for example)
* \todo Filter for intervals (length)
* \todo Memory allocation policy: classic, use a mempool, etc.
*/
 template < class FilteredSimplicialComplexDS
          //, class ArithmeticModifier //only furnishes modifiers and creators
          >
 class Persistent_cohomology {
 public:

  typedef FilteredSimplicialComplexDS               Complex_ds;
  // Data attached to each simplex to interface with a Property Map.
  typedef typename Complex_ds::Simplex_key          Simplex_key;
  typedef typename Complex_ds::Simplex_handle       Simplex_handle;
  typedef typename Complex_ds::Filtration_value     Filtration_value;

  typedef Field_Zp < 2 >                            ArithmeticModifier;
  typedef typename ArithmeticModifier::Element      Arith_element;

// Compressed Annotation Matrix types:
  // Column type
  typedef Cam_column_list < Simplex_key
                          , Arith_element >         Column; // contains 1 set_hook
  // Cell type
  typedef typename Column::Cell                     Cell;   // contains 2 list_hooks
 
  // Remark: constant_time_size must be false because base_hook_cam_h has auto_unlink link_mode
  typedef boost::intrusive::list < Cell
                                 , boost::intrusive::constant_time_size<false> 
                                 , boost::intrusive::base_hook< base_hook_cam_h >    
                                 >                  Hcell;  

  typedef boost::intrusive::set < Column
                                , boost::intrusive::constant_time_size<false> 
                                >                   Cam;
// <------------ try splay_set
 
// Sparse column type for the annotation of the boundary of an element.
  typedef std::vector< std::pair<Simplex_key
                                , Arith_element > > A_ds_type;
// Persistent interval type. The Arith_element field is used for the multi-field framework.
  typedef boost::tuple< Simplex_handle
                      , Simplex_handle
                      , Arith_element >             Persistent_interval;

/** \brief Initializes the Persistent_cohomology class.
  *
  * cpx is a model of FilteredSimplicialComplexDS
  * The Boolean persistence_dim_max is true iff one wants to 
  * compute the persistent homology for the maximal dimension 
  * of faces in the simplicial complex. Default is false.
  */
Persistent_cohomology ( Complex_ds & cpx 
                      , bool         persistence_dim_max = false )
: cpx_    (&cpx)
, dim_max_(cpx.dimension())                   // upper bound on the dimension of the simplices
, ar_pivot_()                                 // initialize the field structure.
, ds_rank_  (cpx_->num_simplices())           // union-find
, ds_parent_(cpx_->num_simplices())           // union-find
, ds_repr_  (cpx_->num_simplices(),NULL)      // union-find -> annotation vectors
, dsets_(&ds_rank_[0],&ds_parent_[0])         // union-find
, cam_()                                      // collection of annotation vectors
, zero_cocycles_()                            // union-find -> Simplex_key of creator for 0-homology
, transverse_idx_()  //, transverse_idx_(cpx_->num_simplices(),NULL) // key -> row
, persistent_pairs_() 
, interval_length_policy(&cpx,0)
, column_pool_(new boost::object_pool< Column > ()) // memory pools for the CAM
, cell_pool_(new boost::object_pool< Cell > ())
{
  if( persistence_dim_max ) { ++dim_max_; }
  // Type Simplex_key must be int
  for(auto sh : cpx_->filtration_simplex_range())
  { dsets_.make_set(cpx_->key(sh)); }
}

// ~Persistent_cohomology()
// { deleting the pools destroys everything they have allocated
// //Clean the remaining columns in the matrix.
//   Column * col_tmp;
//   for(auto cam_it = cam_.begin(); cam_it != cam_.end(); )
//   {
//     col_tmp = &(*cam_it);
//     ++cam_it;
//     delete col_tmp;
//   }
// //Clean the transversal lists
//   for(auto transverse_it = transverse_idx_.begin();
//       transverse_it != transverse_idx_.end(); ++transverse_it)
//   { delete transverse_it->second; }
// delete column_pool_;
// delete cell_pool_;
// }

struct length_interval {
  length_interval ( Complex_ds * cpx
                  , Filtration_value min_length)
  : cpx_(cpx)
  , min_length_(min_length) {}

  bool operator()(Simplex_handle sh1, Simplex_handle sh2)
  { return cpx_->filtration(sh2) - cpx_->filtration(sh1) > min_length_; }

  void set_length(Filtration_value new_length) { min_length_ = new_length; }

  Complex_ds       * cpx_;
  Filtration_value   min_length_;
};
/** \brief Compute the persistent homology of the filtered simplicial
  * complex.
  *
  * Assumes that the filtration provided by the simplicial complex is 
  * valid. Undefined behavior otherwise.*/
void compute_persistent_cohomology ( Filtration_value min_interval_length = 0 )
{
  interval_length_policy.set_length(min_interval_length);

  // Compute all finite intervals
  for( auto sh : cpx_->filtration_simplex_range() )
  {
    int dim_simplex = cpx_->dimension(sh);

    switch(dim_simplex) {
      case 0 :                                              break;
      case 1 : update_cohomology_groups_edge( sh )        ; break;
      default: update_cohomology_groups( sh, dim_simplex ); break;
    }
  }
  // Compute infinite intervals of dimension 0
  Simplex_key key;
  for(auto sh : cpx_->skeleton_simplex_range(0)) //for all vertices
  {
    key = cpx_->key(sh);  
    //if the vertex is representative of its union-find tree and of its
    //connected component as a homology feature, add an interval.
    if( key == dsets_.find_set(key) 
        && zero_cocycles_.find(key) == zero_cocycles_.end() )
    {

      persistent_pairs_.push_back( Persistent_interval ( sh
                                                       , cpx_->null_simplex()
                                                       , ar_pivot_.multiplicative_identity() )
                                  );
    }
  }
  for( auto zero_idx : zero_cocycles_ )
  {
    persistent_pairs_.push_back( Persistent_interval ( cpx_->simplex(zero_idx.second)
                                                     , cpx_->null_simplex()
                                                     , ar_pivot_.multiplicative_identity() )
                                );
  }
// Compute infinite interval of dimension > 0  
// traverse the whole remaining matrix    <----------------- not MF safe
  std::map< Simplex_key,Arith_element > infinite_cocycles;
  typename std::map< Simplex_key,Arith_element >::iterator inf_coc_it;
  for(auto column_it = cam_.begin();
      column_it != cam_.end(); ++column_it)
  {
    for(auto cell_it = column_it->col_.begin();
        cell_it != column_it->col_.end(); ++cell_it)
    {
      inf_coc_it = infinite_cocycles.find(cell_it->key_);
      if( inf_coc_it == infinite_cocycles.end() )
      { infinite_cocycles[cell_it->key_] = 1; }
    }
  }
  //                                      <----------------- not MF safe    inf_coc.second...
  for(auto inf_coc : infinite_cocycles) 
  {
    persistent_pairs_.push_back( Persistent_interval ( cpx_->simplex (inf_coc.first)
                                                     , cpx_->null_simplex()
                                                     , inf_coc.second ) );
  }
}



private:
/** \brief Update the cohomology groups under the insertion of an edge.
  * 
  * The 0-homology is maintained with a simple Union-Find data structure, which
  * explains the existance of a specific function of edge insertions. */
void update_cohomology_groups_edge ( Simplex_handle sigma ) 
{
  Simplex_handle u,v;
  boost::tie(u,v) = cpx_->endpoints(sigma);
  
  Simplex_key ku = dsets_.find_set( cpx_->key(u) ); 
  Simplex_key kv = dsets_.find_set( cpx_->key(v) );

  if(ku != kv ) {        // Destroy a connected component
    dsets_.link(ku,kv);        
    // Keys of the simplices which created the connected components containing
    // respectively u and v. 
    Simplex_key idx_coc_u, idx_coc_v;
    auto map_it_u = zero_cocycles_.find(ku);
    // If the index of the cocycle representing the class is already ku.
    if (map_it_u == zero_cocycles_.end()) { idx_coc_u = ku;               }
    else                                  { idx_coc_u = map_it_u->second; }

    auto map_it_v = zero_cocycles_.find(kv);
    // If the index of the cocycle representing the class is already kv.
    if (map_it_v == zero_cocycles_.end()) { idx_coc_v = kv;               }
    else                                  { idx_coc_v = map_it_v->second; }

    if(idx_coc_u < idx_coc_v) // Kill cocycle [idx_coc_v], which is younger.   
      {
        if(interval_length_policy(cpx_->simplex(idx_coc_v),sigma)) {
          persistent_pairs_.push_back (
              boost::tuple<Simplex_handle,Simplex_handle,Arith_element> ( 
                                                          cpx_->simplex(idx_coc_v)
                                                        , sigma
                                                        , ar_pivot_.multiplicative_identity() 
                                                    )
                                      );
        }
    // Maintain the index of the 0-cocycle alive.
        if( kv != idx_coc_v ) { zero_cocycles_.erase( map_it_v ); }
        if( kv == dsets_.find_set(kv) ) {
          if( ku != idx_coc_u ) { zero_cocycles_.erase( map_it_u ); }
          zero_cocycles_[kv] = idx_coc_u;
        }
      }
    else // Kill cocycle [idx_coc_u], which is younger.
      {
        if(interval_length_policy(cpx_->simplex(idx_coc_u),sigma)) {
          persistent_pairs_.push_back (
              boost::tuple<Simplex_handle,Simplex_handle,Arith_element> ( 
                                                          cpx_->simplex(idx_coc_u)
                                                        , sigma
                                                        , ar_pivot_.multiplicative_identity() 
                                                    )
                                      );
        }
    // Maintain the index of the 0-cocycle alive.
        if( ku != idx_coc_u ) { zero_cocycles_.erase( map_it_u ); }
        if( ku == dsets_.find_set(ku) ) {
          if( kv != idx_coc_v ) { zero_cocycles_.erase( map_it_v ); }
          zero_cocycles_[ku] = idx_coc_v;
        }
      }
    cpx_->assign_key(sigma,cpx_->null_key()); 
  }
  else { // If ku == kv, same connected component: create a 1-cocycle class.
    create_cocycle( sigma, ar_pivot_.multiplicative_identity() ); 
  } 
}


void annotation_of_the_boundary(std::map< Simplex_key, Arith_element > & map_a_ds
                               , Simplex_handle sigma
                               , int dim_sigma )
{
    //traverses the boundary of sigma, keeps track of the annotation vectors,
  // with multiplicity, in a map.
  std::map < Column *, int >                  annotations_in_boundary;
  std::pair < typename std::map< Column *, int >::iterator
            , bool >                          result_insert_bound;
  int sign = 1 - 2 * (dim_sigma % 2); // \in {-1,1} provides the sign in the 
                                      // alternate sum in the boundary.
  Simplex_key key;      Column * curr_col;

  for( auto sh : cpx_->boundary_simplex_range(sigma) )
  {
    key = cpx_->key(sh);
    if( key != cpx_->null_key() ) // A simplex with null_key is a killer, and have null annotation
      {                           // vector.
        // Find its annotation vector
        curr_col = ds_repr_[ dsets_.find_set(key) ];
        if( curr_col != NULL ) 
        { // and insert it in annotations_in_boundary with multyiplicative factor "sign".
          result_insert_bound = 
            annotations_in_boundary.insert(std::pair<Column *,int>(curr_col,sign));  
          if( !(result_insert_bound.second) ) { result_insert_bound.first->second += sign; }
        }
      }
    sign = -sign;
  }
  // Sum the annotations with multiplicity, using a map<key,coeff> 
  // to represent a sparse vector.
  std::pair < typename std::map < Simplex_key, Arith_element >::iterator
            , bool >                                                result_insert_a_ds;

  for( auto ann_ref : annotations_in_boundary ) 
  {    
    if(ann_ref.second != ar_pivot_.additive_identity()) // For all columns in the boundary,
    {                                            
      for( auto cell_ref : ann_ref.first->col_ ) // insert every cell in map_a_ds with multiplicity
      { 
        Arith_element w_y = 
            ar_pivot_.times(cell_ref.coefficient_ , ann_ref.second); //coefficient * multiplicity

        if( w_y != ar_pivot_.additive_identity() ) // if != 0
        {
          result_insert_a_ds = map_a_ds.insert(std::pair< Simplex_key
                                                        , Arith_element >(cell_ref.key_ , w_y));
          if( !(result_insert_a_ds.second) )   //if cell_ref.key_ already a Key in map_a_ds
            { 
              ar_pivot_.plus_equal(result_insert_a_ds.first->second, w_y); 
              if(result_insert_a_ds.first->second == ar_pivot_.additive_identity())
                { map_a_ds.erase(result_insert_a_ds.first); }
            }
        }
      }
    }  
  }

}

/** \brief Update the cohomology groups under the insertion of a simplex.
  * 
  */
void update_cohomology_groups ( Simplex_handle sigma
                              , int dim_sigma )
{
//Compute the annotation of the boundary of sigma:
  std::map< Simplex_key, Arith_element >              map_a_ds;
  annotation_of_the_boundary(map_a_ds, sigma, dim_sigma );

// Update the cohomology groups:
  if( map_a_ds.empty() ) {  // sigma is a creator in all fields represented in ar_pivot_
    if(dim_sigma < dim_max_) { create_cocycle( sigma, ar_pivot_.multiplicative_identity() );}
  }

  else {                    // sigma is a destructor in at least a field in ar_pivot_
 // Convert map_a_ds to a vector
    A_ds_type a_ds; //admits reverse iterators
    for ( auto map_a_ds_ref : map_a_ds )
    { 
      a_ds.push_back( std::pair< Simplex_key
                               , Arith_element> ( map_a_ds_ref.first
                                                , map_a_ds_ref.second ));
    }

   

    Arith_element prod = ar_pivot_.characteristic(); // Product of characteristic of the fields
    for( auto a_ds_rit = a_ds.rbegin(); 
         (a_ds_rit != a_ds.rend()) && (prod != ar_pivot_.multiplicative_identity());
         ++a_ds_rit )
    {
      Arith_element inv_x = ar_pivot_.inverse ( a_ds_rit->second
                                              , prod );     // Modifies "prod"
    /*  
      <----------
      return the charcteristic of fields for which x is invertible.
      put it in the Persistent_interval.
    */
      if( inv_x != ar_pivot_.additive_identity() )
        {
          destroy_cocycle ( sigma
                          , a_ds
                          , a_ds_rit->first
                          , inv_x );
        }
    }
    if( prod != ar_pivot_.multiplicative_identity() && dim_sigma < dim_max_ )
      { create_cocycle( sigma , ar_pivot_.multiplicative_identity(prod) ); }
  }
}

/** \brief Create a new cocycle class.
  *
  * The class is created by the insertion of the simplex sigma.
  * The methods adds a cocycle, representing the new cocycle class,
  * to the matrix representing the cohomology groups.
  * The new cocycle has value 0 on every simplex except on sigma
  * where it worths 1.*/
void create_cocycle ( Simplex_handle sigma
                    , Arith_element x )
{
  Simplex_key key = cpx_->key(sigma);
  // Create a column containing only one cell,
  Column * new_col  = column_pool_->construct(Column(key)); //new Column (key); 
  Cell   * new_cell = cell_pool_->construct(Cell (key, x, new_col)); //new Cell (key, x, new_col);
  new_col->col_.push_back(*new_cell);  
  // and insert it in the matrix, in constant time thanks to the hint cam_.end().
  // Indeed *new_col has the biggest lexicographic value because key is the 
  // biggest key used so far.
  cam_.insert (cam_.end(), *new_col); 
  // Update the disjoint sets data structure.
  Hcell * new_hcell = new Hcell;
  new_hcell->push_back(*new_cell);
  transverse_idx_[key] = new_hcell; //insert the new row
  ds_repr_[key] = new_col;
}

/** \brief Destroy a cocycle class.
*
* The cocycle class is destroyed by the insertion of sigma.
* The methods proceeds to a reduction of the matrix representing 
* the cohomology groups using Gauss pivoting. The reduction zeros-out
* the row containing the cell with highest key in
* a_ds, the annotation of the boundary of simplex sigma. This key
* is "death_key".*/
void destroy_cocycle ( Simplex_handle   sigma
                     , A_ds_type const& a_ds 
                     , Simplex_key      death_key
                     , Arith_element &  inv_x )
{
  // Create a finite persistent interval
  if(interval_length_policy(cpx_->simplex(death_key),sigma)) {
    persistent_pairs_.push_back ( Persistent_interval ( cpx_->simplex(death_key) //creator
                                                      , sigma                    //destructor
                                                      , ar_pivot_.multiplicative_identity() )//fields 
                                );                                   // for which the interval exists 
  }
// <---------- not MF safe



  auto death_key_row = transverse_idx_[death_key]; // Find the beginning of the row.
  std::pair< typename Cam::iterator, bool > result_insert_cam;

  auto row_cell_it = death_key_row->begin();
  while( row_cell_it != death_key_row->end() ) // Traverse all cells in the row at index death_key.
    {
      Arith_element w = ar_pivot_.times_minus( inv_x , row_cell_it->coefficient_ );

      if( w != ar_pivot_.additive_identity() ) 
      { 
        Column * curr_col = row_cell_it->self_col_;         ++row_cell_it;
        // Disconnect the column from the rows in the CAM.
        for( auto col_cell_it = curr_col->col_.begin();
            col_cell_it != curr_col->col_.end(); ++col_cell_it ) 
          { col_cell_it->base_hook_cam_h::unlink(); }
        
        // Remove the column from the CAM before modifying its value
        cam_.erase( cam_.iterator_to(*curr_col) ); 
        // Proceed to the reduction of the column
        plus_equal_column(*curr_col, a_ds, w);

        if( curr_col->col_.empty() ) // If the column is null
        { 
          ds_repr_[ curr_col->class_key_ ] = NULL;  
          column_pool_->free(curr_col); //delete curr_col; 
        }
        else 
        { 
          // Find whether the column obtained is already in the CAM
          result_insert_cam = cam_.insert( *curr_col );
          if ( result_insert_cam.second ) // If it was not in the CAM before: insertion has succeeded
          {
            for ( auto col_cell_it = curr_col->col_.begin();
                  col_cell_it != curr_col->col_.end(); ++col_cell_it ) //re-establish the row links
            { transverse_idx_[ col_cell_it->key_ ]->push_back(*col_cell_it); }
          }
          else // There is already an identical column in the CAM: 
          {    // merge two disjoint sets.
            dsets_.link ( curr_col->class_key_ , 
                          result_insert_cam.first->class_key_ );

            Simplex_key key_tmp = dsets_.find_set( curr_col->class_key_ );
            ds_repr_[ key_tmp ] = &(*(result_insert_cam.first));
            result_insert_cam.first->class_key_ = key_tmp;
            column_pool_->free(curr_col); //delete curr_col;
          }
        }
      }
    else { ++row_cell_it; } // If w == 0, pass.
  }
  // Because it is a killer simplex, set the data of sigma to null_key().
  cpx_->assign_key( sigma, cpx_->null_key() ); 
  transverse_idx_.erase(death_key);
}

/** \brief Assigns target <- target + w * other.*/
void plus_equal_column ( Column & target
                       , A_ds_type const& other //value_type is pair<Simplex_key,Arith_element>
                       , Arith_element w )
{
  auto target_it = target.col_.begin(); auto other_it = other.begin();
  while ( target_it != target.col_.end() && other_it != other.end() )
  {
    if(target_it->key_ < other_it->first) { ++target_it; }
    else {
      if(target_it->key_ > other_it->first) 
      {
        Cell * cell_tmp = cell_pool_->construct(Cell( other_it->first   //key
                                                    , ar_pivot_.additive_identity()
                                                    , &target));
                                   //new Cell( other_it->first   //key
                                   // , ar_pivot_.additive_identity()
                                   // , &target);
        ar_pivot_.plus_times_equal(cell_tmp->coefficient_, other_it->second, w);

        target.col_.insert( target_it, *cell_tmp );

        ++other_it;
      }
      else { //it1->key == it2->key
        //target_it->coefficient_ <- target_it->coefficient_ + other_it->second * w
        ar_pivot_.plus_times_equal( target_it->coefficient_ , other_it->second , w);
        if( target_it->coefficient_ == ar_pivot_.additive_identity() )
        {
          auto tmp_it = target_it;
          ++target_it; ++other_it;   // iterators remain valid
          Cell * tmp_cell_ptr = &(*tmp_it); 
          target.col_.erase(tmp_it); // removed from column
          cell_pool_->free(tmp_cell_ptr);//delete tmp_cell_ptr;       // deleted from memory 
        }
        else { ++target_it; ++other_it; }
      }
    }
  }
}

/** Compare two intervals by length.*/
struct cmp_intervals_by_length {
  cmp_intervals_by_length( Complex_ds * sc ) : sc_ (sc) {}
  bool operator() (  Persistent_interval & p1
                  ,  Persistent_interval & p2 )
  {
    return ( sc_->filtration( get<1>(p1) ) - sc_->filtration( get<0>(p1) ) 
             > sc_->filtration( get<1>(p2) ) - sc_->filtration( get<0>(p2) ) );
  }
  Complex_ds * sc_;
};

public:
/** \brief Outuput the persistence diagram in ostream.
*
* \todo Sort policy, filtered policy.
*/
void output_diagram(std::ostream& ostream = std::cout)
{
//   cmp_intervals_by_length cmp( cpx_ );
//   persistent_pairs_.sort( cmp );
  for(auto pair : persistent_pairs_)
  {
    if(cpx_->filtration(get<0>(pair)) < cpx_->filtration(get<1>(pair))) {
      ostream << cpx_->dimension(get<0>(pair))  << " "
              << cpx_->filtration(get<0>(pair)) << " " 
              << cpx_->filtration(get<1>(pair)) << " " 
              << std::endl;
    }
  }
}

private:
void display_cam() 
{ 
  std::cout << std::endl;
  std::cout << "**** Compressed Annotation Matrix **** \n";

  for(typename Cam::iterator col_it = cam_.begin();
      col_it != cam_.end(); ++col_it ) 
  { col_it->display(); std::cout << std::endl; }

  std::cout << "*************** done ***************** \n";
  std::cout << std::endl;
}



  Complex_ds *         cpx_;
  int                  dim_max_;
  ArithmeticModifier   ar_pivot_;

/** Disjoint sets data structure to link the SimplexDataFilteredSimplicialComplexDS
  * with the compressed annotation matrix.
  * ds_rank_ is a property map Simplex_key -> int, ds_parent_ is a property map 
  * Simplex_key -> simplex_key_t */  
  std::vector< int >                            ds_rank_;  
  std::vector< Simplex_key >                    ds_parent_;
  std::vector< Column * >                       ds_repr_;
  boost::disjoint_sets< int *, Simplex_key * >  dsets_;
/** The compressed annotation matrix fields.*/
  Cam                                           cam_;
/** Dictionary establishing the correspondance between the Simplex_key of
  * the root vertex in the union-find ds and the Simplex_key of the vertex which
  * created the connected component as a 0-dimension homology feature.*/
  std::map<Simplex_key,Simplex_key>             zero_cocycles_;
/** Key -> row. */ 
//  std::vector< Hcell * >                        transverse_idx_;
  std::map< Simplex_key , Hcell * >             transverse_idx_; // PropertyMap Simplex_key -> Hcell *
/** Persistent intervals. */
  std::list< Persistent_interval >              persistent_pairs_; 
  length_interval                               interval_length_policy;


  boost::object_pool< Column > *                column_pool_;
  boost::object_pool< Cell >   *                cell_pool_;
};

#endif // _PERSISTENCECOMPUTATION_SIMPLEXTREE_
