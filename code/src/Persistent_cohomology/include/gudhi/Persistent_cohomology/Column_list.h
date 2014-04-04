/*
 *  Column_list.h
 *  Gudhi
 *
 *  Created by Clément Maria on 15/02/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_COLUMN_LIST_H
#define GUDHI_COLUMN_LIST_H

#include "boost/tuple/tuple.hpp"
#include "boost/intrusive/set.hpp"
#include "boost/intrusive/list.hpp"

template < typename SimplexKey
         , typename ArithmeticElement 
         >
class Cam_column_list;

struct cam_h_tag; // for horizontal traversal in the CAM
struct cam_v_tag; // for vertical traversal in the CAM

typedef boost::intrusive::list_base_hook 
            < boost::intrusive::tag < cam_h_tag > 
            , boost::intrusive::link_mode < boost::intrusive::auto_unlink > //allows .unlink()
            >      base_hook_cam_h;

typedef boost::intrusive::list_base_hook 
            < boost::intrusive::tag < cam_v_tag > 
            , boost::intrusive::link_mode < boost::intrusive::normal_link > //faster hook, less safe
            >      base_hook_cam_v;


template < typename SimplexKey
         , typename ArithmeticElement 
         >
class Cam_matrix_cell 
: public base_hook_cam_h
, public base_hook_cam_v
{
private:
  template < class T1, class T2 > friend class Persistent_cohomology;
  friend class Cam_column_list < SimplexKey , ArithmeticElement >;

  typedef Cam_column_list< SimplexKey, ArithmeticElement > Column;

  Cam_matrix_cell( SimplexKey        key
                 , ArithmeticElement x
                 , Column *          self_col)
  : key_(key)
  , coefficient_(x)
  , self_col_(self_col) {}

  ~Cam_matrix_cell() {}

  SimplexKey                   key_;
  ArithmeticElement            coefficient_;
  Column      *                self_col_;
};





/** \brief Sparse column for the Compressed Annotation Matrix.
*
* The non-zero coefficients of the column are stored in a 
* boost::intrusive::list. Contains a hook to be stored in a
* boost::intrusive::set.
*/
template < typename SimplexKey
         , typename ArithmeticElement >
class Cam_column_list
: public boost::intrusive::set_base_hook 
             < boost::intrusive::link_mode< boost::intrusive::normal_link > > 
{
private:
  template < class T1, class T2 > friend class Persistent_cohomology;

  typedef Cam_matrix_cell < SimplexKey, ArithmeticElement >        Cell;
  typedef boost::intrusive::list < Cell 
                                 , boost::intrusive::constant_time_size<false> 
                                 , boost::intrusive::base_hook< base_hook_cam_v >    
                                 >                                 Col_type;

/** \brief Creates an empty column.*/
  Cam_column_list (SimplexKey key) 
  // : col_()
  // , class_key_(key) 
  {
    class_key_ = key;
    col_ = Col_type();
  }
public:
 /** Copy constructor.*/
 Cam_column_list( Cam_column_list const &other )
 : col_()
 , class_key_(other.class_key_)
 { if(!other.col_.empty()) std::cerr << "Copying a non-empty column.\n"; } 
  private:
  ~Cam_column_list()
  {
  //   typename Col_type::iterator it = col_.begin();
  //   Cell * cell_tmp;
  //   while(it != col_.end())
  //   {
  //     cell_tmp = &(*it);
  //     ++it;
  //     delete cell_tmp; 
  //   }
  }

/** \brief Returns true iff the column is null.*/
  bool is_null() { return col_.empty(); }
/** \brief Returns the key of the representative simplex of
  * the set of simplices having this column as annotation vector
  * in the compressed annotation matrix.*/
  SimplexKey class_key () { return class_key_; }

/** \brief Lexicographic comparison of two columns.*/
friend bool operator< (const Cam_column_list& c1, const Cam_column_list& c2)
  {  
   typename Col_type::const_iterator it1 = c1.col_.begin();
   typename Col_type::const_iterator it2 = c2.col_.begin();
   while(it1 != c1.col_.end() && it2 != c2.col_.end())
   {
      if(it1->key_ == it2->key_)
      { if(it1->coefficient_ == it2->coefficient_) { ++it1; ++it2; }
        else { return it1->coefficient_ < it2->coefficient_; }  }
      else { return it1->key_ < it2->key_; }
   }
  return (it2 != c2.col_.end()); 
  }


  void display()
  {
    for(auto cell : col_) 
    { std::cout << "(" << cell.key_ <<":"<<cell.coefficient_<<") "; }
  }

  Col_type           col_;
  SimplexKey         class_key_;
};

// template <>
// struct lexicographic_compare {
//  typedef Column_cell < SimplexKey, Ring >  Cell;
//   typedef std::list < Cell >                             Col_type;


// bool operator() (Column_list& c1, Column_list& c2)
//   {  
//    typename Col_type::iterator it1 = c1.begin();
//    typename Col_type::iterator it2 = c2.begin();
//    while(it1 != c1.end() && it2 != c2.end())
//    {
//       if(it1->key() == it2->key())
//       { if(it1->coefficient() == it2->coefficient()) { ++it1; ++it2; }
//         else { return it1->coefficient() < it2->coefficient(); }  }
//       else { return it1->key() < it2->key(); }
//    }
//   return (it2 != c2.end()); 
//   }

// };

#endif // GUDHI_COLUMN_LIST_H
