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

//#include "boost/pool.hpp"


/** \brief Computes the persistent cohomology of a filtered simplicial complex.
*
* The implementation is based on the Compressed Annotation Matrix.
*/
 template < class SimplexDataFilteredSimplicialComplexDS
          //, class ArithmeticModifier //only furnishes modifiers and creators
          >
 class Persistent_cohomology {
 public:
// typedefs
// col_idx ?
// typedef boost::intrusive::set< Sparse_col >      boost_set_ann_set; 
// typedef boost_set_ann_set                        ann_set;
// typedef gmp_integer                              ring_elem_t;

  typedef SimplexDataFilteredSimplicialComplexDS    Complex_ds;
/** Data attached to each simplex to interface with a Property Map.*/
  typedef typename Complex_ds::Simplex_key          Simplex_key;
  typedef typename Complex_ds::Simplex_handle       Simplex_handle;
  typedef typename Complex_ds::Filtration_value     Filtration_value;


  typedef Field_Zp < 17 >                      ArithmeticModifier;
  typedef typename ArithmeticModifier::Element Arith_element;

//type of column and matrix
  typedef Cam_column_list < Simplex_key
                          , Arith_element >  Column; // contains 1 set_hook
  typedef typename Column::Cell              Cell;   // contains 2 list_hooks
 
  // Remark: constant_time_size must be false because base_hook_cam_h has
  // auto_unlink link_mode
  typedef boost::intrusive::list < Cell
                                 , boost::intrusive::constant_time_size<false> 
                                 , boost::intrusive::base_hook< base_hook_cam_h >    
                                 >                              Hcell;  

  typedef boost::intrusive::set < Column
                                , boost::intrusive::constant_time_size<false> 
                                >                               Cam;
// try splay_set
 

typedef std::vector< std::pair<Simplex_key, Arith_element > > A_ds_type;


/** \brief Initializes the Persistent_cohomology class.
*
* 
*/
Persistent_cohomology ( Complex_ds & cpx )
: cpx_(&cpx)
, dim_max(cpx.dimension())            // <-- dim_max
//, min_persistence_(1000)
, ar_pivot_()                         //initialize the field structure.
, ds_rank_(cpx_->num_simplices())
, ds_parent_(cpx_->num_simplices())
, ds_repr_(cpx_->num_simplices(),NULL)
, dsets_(&ds_rank_[0],&ds_parent_[0]) //init CAM disjoint sets
, cam_()
, transverse_idx_()
, persistent_pairs_() 
//, column_pool_(new boost::object_pool< Column > ()) // memory pools for the CAM
{
  //valid ? -> transfer in the simplex tree construction?
   //concept for simplex_key_t ?
  for(auto sh : cpx_->filtration_simplex_range())
    {
      dsets_.make_set(cpx_->key(sh));
    }
}

// ~Persistent_cohomology()
// {
// //Clean the remaining columns in the matrix.
//   Column * col_tmp;
//   for(auto cam_it = cam_.begin();
//       cam_it != cam_.end();)
//   {
//     col_tmp = &(*cam_it);
//     ++cam_it;
//     delete col_tmp;
//   }
// //Clean the transversal lists
//   for(auto transverse_it = transverse_idx_.begin();
//       transverse_it != transverse_idx_.end(); ++transverse_it)
//   { delete transverse_it->second; }
// }

/** \todo Parallelize with dimension ?*/
void compute_persistent_cohomology ()
{
  // Show progress: initialization
  boost::progress_display show_progress(cpx_->num_simplices());

  int const dim_max = 10;
  for( auto sh : cpx_->filtration_simplex_range() )
  {
    ++show_progress;

 //   display_cam();
    // if(cpx_->filtration(sh) >= 0.62117) break;

    int dim_simplex = cpx_->dimension(sh);

//    std::cout << "Insert simplex: "; cpx_->display_simplex(sh); std::cout << std::endl;

    switch(dim_simplex) {
      case 0      : break;
      case 1      : update_cohomology_groups_edge( sh );          break;
      case dim_max: update_cohomology_groups_dim_max( sh );       break;
      default     : update_cohomology_groups( sh, dim_simplex );  break; 
    }
  }
}

/** \brief Update the cohomology groups under the insertion of an edge.
* 
* The 0-homology is maintained with a simple Union-Find data structure, which
* explains the existance of a specific function of edge insertions.
*/
void update_cohomology_groups_edge ( Simplex_handle sigma ) 
{
  Simplex_handle u,v;
  boost::tie(u,v) = cpx_->endpoints(sigma);
  
  Simplex_key ku = dsets_.find_set( cpx_->key(u) ); 
  Simplex_key kv = dsets_.find_set( cpx_->key(v) );

  if(ku != kv ) {        // destroys a connected component

    dsets_.link(ku,kv);        
    if(cpx_->filtration(u) < cpx_->filtration(v)) 
      {
         persistent_pairs_.push_back (
            boost::tuple<Simplex_handle,Simplex_handle,Arith_element> ( 
                                                        v
                                                      , sigma
                                                      , ar_pivot_.multiplicative_identity() 
                                                  )
                                    );
      }
    else    
      {
         persistent_pairs_.push_back (
            boost::tuple<Simplex_handle,Simplex_handle,Arith_element> ( 
                                                        u
                                                      , sigma
                                                      , ar_pivot_.multiplicative_identity() 
                                                  )
                                    );
      }
    cpx_->assign_key(sigma,cpx_->null_key()); 
  
  }
  else { // creates a 1-cocycle class
    create_cocycle( sigma, ar_pivot_.multiplicative_identity() ); 
  } 
}

void update_cohomology_groups_dim_max ( Simplex_handle sigma ) 
{
  std::cout << "update_cohomology_groups_dim_max \n";
}

/**
\todo Do we have a problem not considering the orientation?
*/
void update_cohomology_groups ( Simplex_handle sigma
                              , int dim_sigma )
{
   // std::cout << "Enter update_cohomology_groups with      ";
   // cpx_->display_simplex(sigma); std::cout << std::endl;

//Compute the annotation of the boundary of sigma:
  //traverses the boundary of sigma, keeps track of annotations with multiplicities
  // in a map.
//  typedef typename std::map< Column *, int >::iterator ann_in_bound_iterator;

  std::map < Column *, int >                       annotations_in_boundary;
  std::pair < typename std::map< Column *, int >::iterator
            , bool >                               result_insert_bound;
  
  int sign = 1 - 2 * (dim_sigma % 2); // \in {-1,1} provides the sign for the 
                                      // alternate sum in the boundary
  Simplex_key key;  
  Column * curr_col;
  
  for( auto sh : cpx_->boundary_simplex_range(sigma) )
  {
 //    std::cout << "      --- "; cpx_->display_simplex(sh); std::cout << "     ";

    key = cpx_->key(sh);

//    std::cout << " (" << key << ") ";

    if( key != cpx_->null_key() ) // having a null_key means be a killer simplex, and killer
      {                           // simplices have null annotation.
        //find its annotation vector
        curr_col = ds_repr_[ dsets_.find_set(key) ];
        if( curr_col != NULL ) 
        { //and insert it in annotations_in_boundary with coefficient sign
          result_insert_bound = 
            annotations_in_boundary.insert(std::pair<Column *,int>(curr_col,sign));  
          if( !(result_insert_bound.second) ) { result_insert_bound.first->second += sign; }
     
    //      curr_col->display();

        }

    //    else { std::cout << " curr_col = NULL ";}

      }


    if( sign == 1 ) sign = -1;
    else            sign =  1;
  

//    std::cout << std::endl;
  }


// std::cout << "   Map of boundary annotations: \n";
//   for( auto ann_ref : annotations_in_boundary ) 
// {
//   std::cout << "      " << ann_ref.second << "   * ";
//   ann_ref.first->display();
//   std::cout << std::endl;
// }

  // sums the annotations with multiplicity, using a map<key,coeff> 
  // to represent a sparse vector.
  std::map< Simplex_key, Arith_element >                       map_a_ds;
  std::pair < typename std::map < Simplex_key
                                , Arith_element >::iterator
            , bool >                                           result_insert_a_ds;

  for( auto ann_ref : annotations_in_boundary ) 
  {    
    if(ann_ref.second != ar_pivot_.additive_identity())
    {                                            //for all columns in the boundary
      for( auto cell_ref : ann_ref.first->col_ ) //insert every cell in map_a_ds with multiplicity
      { 
        Arith_element w_y = 
            ar_pivot_.times(cell_ref.coefficient_ , ann_ref.second); //coefficient * multiplicity

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

// Updates the cohomology groups:
  if( map_a_ds.empty() )  // sigma is a creator in all fields represented in ar_pivot_
   { create_cocycle( sigma, ar_pivot_.multiplicative_identity() ); }

  else {                  // sigma is a destructor in at least a field in ar_pivot_
 // Convert map_a_ds to a vector
    A_ds_type a_ds; //admits reverse iterators
    for ( auto map_a_ds_ref : map_a_ds )
    { 
      a_ds.push_back( std::pair< Simplex_key
                               , Arith_element> ( map_a_ds_ref.first
                                                , map_a_ds_ref.second ));
    }




    // std::cout << "    Value of A_ds: ";
    // for(auto a_ds_p : a_ds)
    //   { std::cout <<"("<<a_ds_p.first<<":"<<a_ds_p.second<<") ";}
    // std::cout<<std::endl;




    Arith_element prod = ar_pivot_.characteristic(); // product of characteristic of the fields

    for( auto a_ds_rit = a_ds.rbegin(); 
        (a_ds_rit != a_ds.rend()) && (prod != ar_pivot_.multiplicative_identity());
        ++a_ds_rit )
    {
      Arith_element inv_x = ar_pivot_.inverse ( a_ds_rit->second
                                              , prod );        // <- Careful : modifies prod
      
      //std::cout<<"   inverse of lowest non-zero coeff in a_ds = " << inv_x << std::endl;

      if( inv_x != ar_pivot_.additive_identity() )
        {
          destroy_cocycle ( sigma
                          , a_ds
                          , a_ds_rit->first
                          , inv_x ); 
        }
    }

    if( prod != ar_pivot_.multiplicative_identity() )
      { create_cocycle( sigma , ar_pivot_.multiplicative_identity(prod) ); }
  }
}

/** \brief Creates a new cocycle class.
  *
  * The class is created by the insertion of the simplex sigma.
  * The methods adds a cocycle, representing the new cocycle class,
  * to the matrix representing the cohomology groups.
  * The new cocycle has value 0 on every simplex except on sigma
  * where it worths 1.*/
void create_cocycle ( Simplex_handle sigma
                    , Arith_element x )
{

   // std::cout << "Creator of cocycle: " << cpx_->key(sigma);
   // std::cout << std::endl;


  Simplex_key key = cpx_->key(sigma);
  //creates a column containing only one cell.
  Column * new_col  = new Column (key); 
  Cell   * new_cell = new Cell (key, x, new_col);
  new_col->col_.push_front(*new_cell);  


  //insert it in the matrix, in constant time thanks to the hint cam_.end().
  //indeed *new_col has the biggest lexicographic value because key is the 
  //biggest key used so far.
  cam_.insert (cam_.end(), *new_col); 

  //update the disjoint sets data structure.
  Hcell * new_hcell = new Hcell();
  new_hcell->push_front(*new_cell);
  transverse_idx_[key] = new_hcell; //insert the new row

  ds_repr_[key] = new_col;
}

/** \brief Destroys a cocycle class.
*
* The cocycle class is destroyed by the insertion of sigma.
* The methods proceeds to a reduction of the matrix representing 
* the cohomology groups using Gauss pivoting. The reduction zeros-out
* the row containing lowest_cell, wioth is the cell with highest key in
* a_ds, the annotation of the boundary of simplex sigma.*/
void destroy_cocycle ( Simplex_handle sigma
                     , A_ds_type const& a_ds 
                     , Simplex_key death_key
                     , Arith_element & inv_x )
{

   // std::cout << "Destructor of cocycle: " << death_key << "\n";

  //set the data of sigma to null_key()
  cpx_->assign_key( sigma, cpx_->null_key() ); 
  //create a persistent pair
  persistent_pairs_.push_back (
      boost::tuple<Simplex_handle,Simplex_handle,Arith_element> ( 
                                      cpx_->simplex(death_key)      //creator
                                    , sigma                                //destructor
                                    , ar_pivot_.multiplicative_identity() )//fields coeff in for
                              );                                           //which the interval exists 

  auto death_key_row = transverse_idx_.find( death_key ); //find the beginning of the row
  std::pair< typename Cam::iterator, bool > result_insert_cam;

  //in order to access the Hcell containg a given cell
  //typedef boost::intrusive::circular_list_algorithms< base_hook_cam_h > algoHcell;

  auto row_cell_it = death_key_row->second->begin();

  while( row_cell_it != death_key_row->second->end() ) // traverse all cells in the row at index 
    {                                                  // the index of lowest_cell
      Arith_element w = ar_pivot_.times_minus( inv_x , row_cell_it->coefficient_ );


  // std::cout << "D \n";

  //   std::cout << "----------- "<< row_cell_it->coefficient_ << "      w = " << w << std::endl;

  // std::cout << "E \n";



      if( w != ar_pivot_.additive_identity() ) 
      { //disconnect the column from the CAM
        Column * curr_col = row_cell_it->self_col_;  ++row_cell_it;


  // std::cout << "F \n";



        for( auto col_cell_ref : curr_col->col_ ) 
          {
//           algoHcell::unlink( &col_cell_ref ); 
            col_cell_ref.base_hook_cam_h::unlink(); 
          }

  // std::cout << "G \n";


        //remove the column from the CAM before modifying its value
        cam_.erase( cam_.iterator_to(*curr_col) ); 
      
  // std::cout << "H \n";


        //proceed to the reduction of the column

      // std::cout << "Plus equal: ";
      // curr_col->display();
      // std::cout << "      ---> "; 
     
    //  curr_col->display(); std::cout << "  +  " << w << " * ";
    // for(auto a_ds_p : a_ds)
    //   { std::cout <<"("<<a_ds_p.first<<":"<<a_ds_p.second<<") ";}
    // std::cout<< "  =  " << std::endl << "            ";


        plus_equal_column(*curr_col, a_ds, w);   // <- w is an Arith_element, not an int
     

       // curr_col->display();
       // std::cout << std::endl;


  // std::cout << "I \n";





      //find whether the column obtained is already in the CAM
      if( curr_col->col_.empty() ) // if the column is null
        { 


  //           std::cout << "   !!!!!!!!!!!!!!!! Empty !!!!!!!!!!!!!!!! \n";

          if(curr_col->class_key_ != dsets_.find_set(curr_col->class_key_))
            { std::cerr << " INVALID COLUMN CLASS KEY \n \n \n ";}

          ds_repr_[ curr_col->class_key_ ] = NULL;  
          delete curr_col; 
        }
      else 
        { 

             // std::cout << "Not empty \n" << "                  ";
             // curr_col->display();
             // std::cout << std::endl;

          result_insert_cam = cam_.insert( *curr_col );
          if ( result_insert_cam.second ) // was not in the CAM before: insertion has succeeded
          {
            for ( auto col_cell_ref : curr_col->col_ ) //re-establish the row links
            { transverse_idx_[ col_cell_ref.key_ ]->push_front(col_cell_ref); }
          }
          else
            { // already in the CAM, merges two disjoint sets 
             
            //  std::cout << "       already in CAM \n";

            // std::cout << "       curr_col key = " << curr_col->class_key_ << "   and other = " 
            //           << result_insert_cam.first->class_key_ << std::endl;

              dsets_.link ( curr_col->class_key_ , 
                            result_insert_cam.first->class_key_ );

              // std::cout << "   Result link is " << dsets_.find_set( curr_col->class_key_ )
              //           << "     and    " << dsets_.find_set( result_insert_cam.first->class_key_ )
              //           << std::endl;

              Simplex_key key_tmp = dsets_.find_set( curr_col->class_key_ );

              ds_repr_[ key_tmp ] = &(*(result_insert_cam.first));
              
              result_insert_cam.first->class_key_ = key_tmp;

              delete curr_col;
            }
        }
      }
    else { ++row_cell_it; } // w == 0
  }
}





/** \brief Assigns target <- target + w * other. */
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
        Cell * cell_tmp = new Cell( other_it->first   //key
                                   , ar_pivot_.additive_identity()
                                   , &target);
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
          delete tmp_cell_ptr;       // deleted from memory 
        }
        else { ++target_it; ++other_it; }
      }
    }
  }
}








/** Compare two intervals by length.*/
struct cmp_intervals_by_length {

  cmp_intervals_by_length( Complex_ds * sc ) : sc_ (sc) {}

  bool operator() (  boost::tuple< Simplex_handle, Simplex_handle, Arith_element > & p1
                  ,  boost::tuple< Simplex_handle, Simplex_handle, Arith_element > & p2 )
  {
    return ( sc_->filtration( get<1>(p1) ) - sc_->filtration( get<0>(p1) ) 
             > sc_->filtration( get<1>(p2) ) - sc_->filtration( get<0>(p2) ) );
  }

  Complex_ds * sc_;
};

void output_diagram(std::ostream& ostream = std::cout)
{
//  std::cout << "enter output_diagram \n";
  std::cout << "Number of pairs = " << persistent_pairs_.size() << std::endl;
  cmp_intervals_by_length cmp( cpx_ );
  persistent_pairs_.sort( cmp );
  for(auto pair : persistent_pairs_)
  {
    // ostream << cpx_->filtration(get<0>(pair)) << " " 
    //         << cpx_->filtration(get<1>(pair)) << "    " 
    //         << get<2>(pair) << std::endl;

    cpx_->display_simplex(get<0>(pair));
    std::cout << "   -   ";
    cpx_->display_simplex(get<1>(pair));
    std::cout << "       in " << get<2>(pair) << std::endl;
  }
}




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
  int                  dim_max;
  //Filtration_value     min_persistence_; <- use a predicate instead?
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
  std::map< Simplex_key , Hcell * >             transverse_idx_;

//  boost::object_pool< Column > column_pool_;
//  boost::object_pool< Cell >   cell_pool_;


//more compact when using a single field?
  std::list< boost::tuple<Simplex_handle,Simplex_handle,Arith_element> > persistent_pairs_;

};




// /** Compare two intervals by length.*/
// template < class SimplexDataSimplicialComplexDS
//          , typename Arith_element >
// struct cmp_intervals_by_length {
//   typedef typename SimplexDataSimplicialComplexDS::Simplex_handle Simplex_handle;

//   cmp_intervals_by_length( SimplexDataSimplicialComplexDS * sc ) : sc_ (sc) {}

//   bool operator() (  boost::tuple< Simplex_handle, Simplex_handle, Arith_element > & p1
//                   ,  boost::tuple< Simplex_handle, Simplex_handle, Arith_element > & p2 )
//   {
//     return ( sc_->filtration( get<1>(p1) ) - sc_->filtration( get<0>(p1) ) 
//              > sc_->filtration( get<1>(p2) ) - sc_->filtration( get<0>(p2) ) );
//   }

//   SimplexDataSimplicialComplexDS * sc_;
// };

// /** Output the diagram as a set of intervals. 
//   * Sort the intervals by length first.*/
// template < class Sc >
// std::ostream& operator<< ( std::ostream& ostream
//                          , Persistent_cohomology< Sc > & pcoh)
// {
//   cmp_intervals_by_length<Sc, typename Persistent_cohomology< Sc >::Arith_element > cmp(pcoh.cpx_);
//   pcoh.persistent_pairs_.sort( cmp );
//   for(auto pair : pcoh.persistent_pairs_)
//   {
//     std::cout << get<0>(pair) << " " << get<1>(pair) << "    " 
//               << get<2>(pair) << std::endl;
//   }
// }













#endif // _PERSISTENCECOMPUTATION_SIMPLEXTREE_
