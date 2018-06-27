#include <iostream>

/** Represents an edge for zigzag filtrations for flag complexes. 
  * The edge must have two endpoints, encoded by Vertex_handles, a filtration 
  * value and a type (insertion or deletion) represented by a bool.
  *
  * A sequence of such edges represents a full flag zigzag filtration.
  *
  * An edge with u_ == v_ represents a vertex.
  */
template< class FilteredComplex >
class Zigzag_edge {
public:
  Zigzag_edge( typename FilteredComplex::Vertex_handle u
             , typename FilteredComplex::Vertex_handle v
             , typename FilteredComplex::Filtration_value fil
             , bool type)
  : u_(u), v_(v), fil_(fil), type_(type) {}

/* Returns vertex with smaller label. */
  typename FilteredComplex::Vertex_handle    u()    { return u_; }
/* Returns vertex with bigger label. */
  typename FilteredComplex::Vertex_handle    v()    { return v_; }
/* Returns the filtration value of the edge. */
  typename FilteredComplex::Filtration_value fil()  { return fil_; }
/* Returns true if insertion of the edge, false if removal. */
  bool                                       type() { return type_; }

  bool operator==(const Zigzag_edge &e) const {
    return ( (e.u_ == u_) && (e.v_ == v_) && 
             (e.fil_ == fil_) && (e.type_ == type_) );
  }

private:
  typename FilteredComplex::Vertex_handle    u_;
  typename FilteredComplex::Vertex_handle    v_;
  typename FilteredComplex::Filtration_value fil_;
  bool                                       type_;
};

/** 
  * Iterator over a flag zigzag filtration implicitely 
  * represented by a list of Zigzag_edges. 
  *
  * Given an empty FlagZigzagFiltrationComplex and a range of insertions and 
  * deletion of edges, the iterator add/remove on the fly the range of edges and 
  * expand the complex in consequence. It traverses all the newly added/removed 
  * simplices (induced by the new edge) before doing further modifications.
  *
  */
template< class FlagZigzagFilteredComplex >
class Flagzigzag_simplex_iterator 
: public boost::iterator_facade<
            Flagzigzag_simplex_iterator<FlagZigzagFilteredComplex>
          , typename FlagZigzagFilteredComplex::Simplex_handle
          , boost::forward_traversal_tag >
{
  public:
    typedef typename FlagZigzagFilteredComplex::Simplex_handle Simplex_handle;
    typedef typename FlagZigzagFilteredComplex::Edge_type      Edge_type;

    Flagzigzag_simplex_iterator() //any end iterator
    : cpx_(NULL)
    // , zigzag_edge_filtration_(NULL)
    , counter_insert(0) 
    {}

//cpx must be empty
    Flagzigzag_simplex_iterator( FlagZigzagFilteredComplex * cpx 
                               , std::vector< Edge_type >  * zz_edge_fil_ptr
                               , int                         dim_max )
    {
      zigzag_edge_filtration_ = zz_edge_fil_ptr;
      dim_max_                = dim_max;
      are_we_done             = false;
      cpx_                    = cpx;
      counter_insert          = 0;
      partial_zzfil_          = std::vector< Simplex_handle >(); //TODO?
      edge_it_                = zigzag_edge_filtration_->begin();
      if(edge_it_ == zigzag_edge_filtration_->end()) 
      { cpx_ = NULL; return; } //end() iterator
      
      //add the first edge
      arrow_direction_ = edge_it_->type(); //must be true, i.e., and insertion
      cpx_->flag_add_edge( edge_it_->u(), edge_it_->v(), edge_it_->fil()
                                , dim_max_, partial_zzfil_);
      sh_it_ = partial_zzfil_.begin();
      ++edge_it_;
      for(auto & sh : partial_zzfil_) 
      { sh->second.assign_key(counter_insert); ++counter_insert; } 
    }

//User-defined copy constructor
    Flagzigzag_simplex_iterator(const Flagzigzag_simplex_iterator& other )
    : cpx_(other.cpx_)
    , zigzag_edge_filtration_(other.zigzag_edge_filtration_)
    , dim_max_(other.dim_max_)
    , partial_zzfil_(other.partial_zzfil_)
    , sh_it_(partial_zzfil_.begin())
    , edge_it_(other.edge_it_)
    , arrow_direction_(other.arrow_direction_)
    , counter_insert(other.counter_insert)
    , are_we_done(other.are_we_done) {}

    bool arrow_direction() { return arrow_direction_; }

  private:
    friend class boost::iterator_core_access;

    bool equal(Flagzigzag_simplex_iterator const& other) const {
      if(cpx_ == NULL) { return (other.cpx_ == NULL); }      
      return ( cpx_     == other.cpx_     && 
               edge_it_ == other.edge_it_ &&
               sh_it_   == other.sh_it_ );
    }

    Simplex_handle & dereference() const {
      return *sh_it_;
    }

    void increment() 
    {
      ++sh_it_;
      if(sh_it_ == partial_zzfil_.end()) //add or remove the next edge
      { //check if we have reached the end of a sequence of backward arrows, 
        //associated to the removal of an edge. If so, we remove effectively 
        //the simplices from the complex.
        if(!arrow_direction_) //need to effectively remove the simplices we have just considered.
        { 
          //effectively remove all simplices from partial_zzfil_; must be sorted 
          cpx_->remove_maximal_simplices(partial_zzfil_);

          //The simplices in partial_zzfil_ come by decreasing keys, hence 
          //are all maximal when removing from left to right. 
          //FlagZigzagFilteredComplex::Dictionary must not invalidate iterators
          //when inserting and removing (e.g., std::map<,>).
          //IF we want to maintain the validity of Simplex_handle (i.e. map 
          //iterators) during removals, even when using boost::flat_map. To do 
          //so, we add a sorting procedure, maintaining the maximality 
          //property, so as we remove Nodes in a Siblings::members() from 
          //right to left.
          // sort( partial_zzfil_.begin(), partial_zzfil_.end()
          //     , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
          //       if(sh1->first != sh2->first) {return sh1->first > sh2->first;}
          //       return sh1->second.key() > sh2->second.key();
          //     });

        //   for( auto sh_it = partial_zzfil_.begin();
        //             sh_it != partial_zzfil_.end(); ++sh_it) 
        //   { 
        //     // (*sh_it)->second.unlink_hooks();
        //     // (*sh_it)->second.assign_key(-21);
        //     cpx_->remove_maximal_simplex(*sh_it); //modify the complex 
        //   } 
        // }


        }


        partial_zzfil_.clear();
     
        //if all edges have been considered:
        if(edge_it_ == zigzag_edge_filtration_->end()) 
        { 
          if(are_we_done) { cpx_ = NULL; return; } //set iterator to end() position 
          else {//no edge left, consider simplices remaining in the complex 
            are_we_done = true;//happens once
            //fills up zz_partial with the remaining simplices in complex
            cpx_->flag_lazy_empty_complex(partial_zzfil_); 
            arrow_direction_ = false; //only backward arrows now

            sort( partial_zzfil_.begin(), partial_zzfil_.end()
                , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
                    return sh1->second.key() > sh2->second.key();
                });

            sh_it_ = partial_zzfil_.begin();
            return;
          }
        }
        //partial_zzfil_ is empty
        if( edge_it_->type() ) { //forward arrow //modify the complex
          cpx_->flag_add_edge( edge_it_->u(), edge_it_->v()
                             , edge_it_->fil()
                             , dim_max_, partial_zzfil_ );
          arrow_direction_ = true; //the arrow is forward

          //flag_add_edge output a SORTED sequence of simplices
          for(auto & sh : partial_zzfil_) //set key values
          { sh->second.assign_key(counter_insert); ++counter_insert; }
        }
        else { //backward arrow
          cpx_->flag_lazy_remove_edge( edge_it_->u(), edge_it_->v()
                                     , partial_zzfil_ ); //does not modify cpx
          arrow_direction_ = false; //the arrow is backward

          sort( partial_zzfil_.begin(), partial_zzfil_.end()
              , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
                  return sh1->second.key() > sh2->second.key();
              });
        }
       //partial_zzfil_ contains at least the new edge
        sh_it_ = partial_zzfil_.begin(); 
        ++edge_it_;
      }
    }
  


  //complex getting modified
  FlagZigzagFilteredComplex                                            * cpx_; 
/* List of insertion and deletion of edges representing the flag zigzag fil.*/
  std::vector< Edge_type >                           * zigzag_edge_filtration_;
/* Maximal dimension of the flag complex. */
  int                                                                 dim_max_;
/* part of the zz filtration constructed by the last edge insertion. 
 * When reaching the end of it, clear it, insert a new edge via the edge_it_++ 
 * and compute a new bunch of simplices of the zz filtration. */
  typename std::vector< Simplex_handle >           partial_zzfil_;
  //current simplex in partial_zzfil_
  typename std::vector< Simplex_handle >::iterator sh_it_;
  //iterator in the range of edges; points to the next edge to insert or remove
  typename std::vector< Edge_type >::iterator      edge_it_;
  //true if the simplices in partial_zzfil_ are insertions, and false if deletions
  bool                                             arrow_direction_;
  //counts the total number of insertions in the zigzag, used to assign keys
  int                                              counter_insert;
  //true iff we are finishing emptying the complex
  bool                                             are_we_done;
};
