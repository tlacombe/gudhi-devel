// #include <iostream>

/** Represents an edge for zigzag filtrations for flag complexes. 
  * The edge must have two endpoints, encoded by Vertex_handles, a filtration 
  * value and a type (insertion or deletion).
  *
  * A sequence of such edges represents a full flag zigzag filtration.
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

private:
  typename FilteredComplex::Vertex_handle    u_;
  typename FilteredComplex::Vertex_handle    v_;
  typename FilteredComplex::Filtration_value fil_;
  bool                                       type_;
};

/** Iterator over a flag zigzag filtration implicitely 
  * represented by a list of Zigzag_edges. 
  */
template< class FlagZigzagFilteredComplex >
class Flagzigzagfiltration_simplex_iterator 
: public boost::iterator_facade<
            Flagzigzagfiltration_simplex_iterator<FlagZigzagFilteredComplex>
          , typename FlagZigzagFilteredComplex::Simplex_handle
          , boost::forward_traversal_tag >
{
  public:
    typedef typename FlagZigzagFilteredComplex::Simplex_handle Simplex_handle;
    typedef typename FlagZigzagFilteredComplex::Edge_type      Edge_type;

    Flagzigzagfiltration_simplex_iterator() //any end iterator
    : cpx_(NULL)
    , counter_insert(0) {}

    Flagzigzagfiltration_simplex_iterator(FlagZigzagFilteredComplex * cpx)
    {
      are_we_done    = false;
      cpx_           = cpx;
      counter_insert = 0;
      partial_zzfil_ = std::vector< Simplex_handle >(); //TODO?
      edge_it_       = cpx->zigzag_edge_filtration_->begin();
      if(edge_it_ == cpx->zigzag_edge_filtration_->end()) 
      { cpx_ = NULL; return; } //end() iterator
      //add the first edge
      arrow_direction_ = edge_it_->type();

      cpx_->flag_zigzag_add_edge( edge_it_->u(), edge_it_->v(), edge_it_->fil()
                                , cpx_->dim_max_, partial_zzfil_);

      sh_it_ = partial_zzfil_.begin();
      ++edge_it_;

      for(auto & sh : partial_zzfil_) 
      { sh->second.assign_key(counter_insert); ++counter_insert; } 
    }

    bool arrow_direction() { return arrow_direction_; }

  private:
    friend class boost::iterator_core_access;

    bool equal(Flagzigzagfiltration_simplex_iterator const& other) const {
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
      if(sh_it_ == partial_zzfil_.end())
      {
        //check if we have reached the end of a sequence of backward arrows, 
        //associated to the removal of an edge. If so, we remove effectively 
        //the simplices from the complex.
        if(!arrow_direction_) 
        { //each simplex is maximal in the simplex tree due to filtration order
          
          //reorder partial_zzfil_ so as iterators do not get invalidated when suppressing Simplex_handles.
          sort( partial_zzfil_.begin(), partial_zzfil_.end()
              , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
                //we want the rightmost elements of the flat_maps to be removed first, so as the Simplex_handle on the left, that may also be in partial_zzfil_, do not get invalidated.
                  if(sh1->first != sh2->first) {return sh1->first > sh2->first;}
                  return sh1->second.key() > sh2->second.key();
              });

          for( auto sh_it = partial_zzfil_.begin();
               sh_it != partial_zzfil_.end(); ++sh_it) 
          { 
            (*sh_it)->second.unlink_hooks();
            (*sh_it)->second.assign_key(-21);
            cpx_->remove_maximal_simplex(*sh_it); //modify the complex 
          } 
        }
        partial_zzfil_.clear();
        //if all edges have been considered:
        if(edge_it_ == cpx_->zigzag_edge_filtration_->end()) 
        { 
          if(are_we_done) 
          { cpx_ = NULL; return; } //set iterator to end() position 
          else //no edge left, but there may be simplices remaining in the complex, like vertices. 
          {
            are_we_done = true;
            //fills up zz_partial with the remaining simplices in complex
            cpx_->flag_zigzag_lazy_empty_complex(partial_zzfil_); 
            arrow_direction_ = false; //only backward arrows now

            sort( partial_zzfil_.begin(), partial_zzfil_.end()
                , [](Simplex_handle sh1, Simplex_handle sh2)->bool {
                    return sh1->second.key() > sh2->second.key();
                });

            sh_it_ = partial_zzfil_.begin();
            return;
          }
        }

        if( edge_it_->type() ) //forward arrow
        { //modify the complex
          cpx_->flag_zigzag_add_edge( edge_it_->u(), edge_it_->v()
                                    , edge_it_->fil()
                                    , cpx_->dim_max_, partial_zzfil_ );
          arrow_direction_ = true; //the arrow is forward

          for(auto & sh : partial_zzfil_) //set key values
          { sh->second.assign_key(counter_insert); ++counter_insert; }
        }
        else { //backward arrow
          cpx_->flag_zigzag_lazy_remove_edge( edge_it_->u(), edge_it_->v()
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

//complex getting manipulated
  FlagZigzagFilteredComplex                      * cpx_; 
  //part of the zz filtration constructed by the last edge insertion. When reaching
  //the end of it, clear it, insert a new edge via the edge_it_++ and compute a new 
  //bunch of simplices of the zz filtration.
  typename std::vector< Simplex_handle >           partial_zzfil_;
  //current simplex in current_zz_fil_
  typename std::vector< Simplex_handle >::iterator sh_it_;
  //iterator in the set of edges stored in cpx_. points to the next edge to insert/remove
  typename std::vector< Edge_type >::iterator      edge_it_;
  //true if the simplices in partial_zzfil_ are insertions, and false if deletions
  bool                                             arrow_direction_;
  //counts the total number of insertions in the zigzag
  int                                              counter_insert;
  //true iff we are finishing emptying the complex
  bool                                             are_we_done;
};
