
/*---------------------------------------------------------------------------*/

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
    if(sib_->oncles() != NULL)    {    sib_ = sib_->oncles();}
    else                        {    sib_ = NULL    ;}    
    return *this;
  }

 private:
  Siblings *    sib_;
  Vertex            v_;
};
/*------------------------------------*/
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

/*---------------------------------------------------------------------------*/

template<class N>
class Flag_simplex_tree<N>::Boundary_simplex_iterator {
 public:
 Boundary_simplex_iterator() :
  sib_(NULL)
    {}

  /**
   */
 Boundary_simplex_iterator(    Flag_simplex_tree * st,
			       Simplex_handle sh) :
  suffix_(),
    st_(st)
    {
      if(st == NULL) //end()
	{    sib_ = NULL;    }
      else 
	{
	  last_            = sh->first;
	  Siblings * sib    = sh->second.self_siblings(last_);
	  next_            = sib->parent();
	  sib_            = sib->oncles();       /** \todo check if NULL*/
	  sh_                = sib_->find(next_);
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
  Vertex                            last_   ; //last vertex of the simplex
  Vertex                            next_   ; //next vertex to push in suffix_
  std::vector< Vertex >             suffix_ ;
  Siblings                    *     sib_    ; //where the next search will start from
  Simplex_handle                    sh_     ; //current Simplex_handle in the boundary
  Flag_simplex_tree           *     st_     ; //simplex containing the simplicial complex
};
/*-----------------------------------*/
/**
 * Range over the simplices in a boundary.
 * .begin() and .end() return Boundary_simplex_iterator
 * type.
 */
template<class N>
class Flag_simplex_tree<N>::Boundary_simplex_range {
 public:
 Boundary_simplex_range(Flag_simplex_tree *    st,
			Simplex_handle                sh) :
  st_(st),
    sh_(sh)
    {}

  Boundary_simplex_iterator begin()
  {return Boundary_simplex_iterator(st_,sh_);    }

  Boundary_simplex_iterator end()
  {return Boundary_simplex_iterator(NULL,sh_);}

 private:
  Flag_simplex_tree        *     st_;
  Simplex_handle                sh_;
};

/*---------------------------------------------------------------------------*/

template<class N>
class Flag_simplex_tree<N>::Filtration_simplex_range {
 public:
 Filtration_simplex_range(Flag_simplex_tree * st) :
  st_(st)
  {}
  
  Filtration_simplex_iterator begin()
  { return st_->filtration_vector().begin();}

  Filtration_simplex_iterator end()
  {st_->filtration_vector().end();}

 private:
  Flag_simplex_tree      *    st_;
};

/*---------------------------------------------------------------------------*/

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
  Simplex_handle        sh_;
  Flag_simplex_tree  *  st_;
};
/*-----------------------------------*/
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

/*---------------------------------------------------------------------------*/
