
template<class N >
Flag_simplex_tree<N>::
Flag_simplex_tree() :
  gt_(new N()),
  rho_max_(0),
  nb_vertices_(0),
  size_cpx_(0),
  root_(),
  filtration_vect_()
{}

template<class N>
void 
Flag_simplex_tree<N>::print(Simplex_handle sh,
          std::ostream& os)
{
  Simplex_vertex_range svr = simplex_vertex_range(sh);
  for(Simplex_vertex_iterator it = svr.begin();
      it != svr.end(); ++it)
    {os << *it << " ";}
  os << std::endl;
}

template<class N >
typename Flag_simplex_tree<N>::Simplex_handle 
Flag_simplex_tree<N>::find(std::vector < Vertex > & s)
{
  if(s.size() == 0) std::cerr << "Empty simplex \n";
  if(s.size() == 1) std::cerr << "Vertex \n";

  //for(typename std::vector<Vertex>::iterator it_ = s.begin(); 
  //it_ != s.end();++it_)
  //{std::cout << *it_ << " ";}
  //std::cout << " \n";

  typename std::vector< Vertex >::iterator it = s.begin();
  if(! root_[ *it ].has_children( *it )) { }    // not there
  Siblings * for_sib = root_[ *it ].children();
  ++it;
  Simplex_handle sh = for_sib->find( *it ); //some stop condition here
  ++it;
  for( ; it != s.end(); ++it)
    {
      if(! sh->second.has_children(sh->first)) std::cerr << "Not here \n";    //must create a st.end() Simplex_handle...
      sh = sh->second.children()->find(*it);// return some false if not here...
    }
  return sh;
}

template<class N>
void 
Flag_simplex_tree< N >::init(//Point_range_sc&point_range,
           int dim_max,
           Filtration_value rho_max)
{
  rho_max_ = rho_max;
  nb_vertices_ = gt_->nb_elements();
  // Insert all edges
  root_ = std::vector< Node >( gt_->nb_elements(), Node() );
  for(Vertex_iterator v_it = gt_->vertex_range().begin();
      v_it != gt_->vertex_range().end(); ++v_it)
    {
      Neighbor_vertex_range n_range(gt_,*v_it,rho_max);
      for(Neighbor_vertex_iterator n_it = n_range.begin();
    n_it != n_range.end(); ++n_it)
  {
    if(*v_it < *n_it) 
      {
        if(! root_[*v_it].has_children(*v_it)) 
    {root_[*v_it].assign_children(new Siblings(NULL,*v_it)); }
        root_[*v_it].children()->insert(*n_it,gt_->distance(*v_it,*n_it));
      }
  }
    }
  // Update size of the complex
  size_cpx_ += root_.size();
  int v = 0;
  for(std::vector< Node >::iterator r_it = root_.begin();
      r_it != root_.end(); ++r_it,++v)
    { if(r_it->has_children(v)) {size_cpx_ += r_it->children()->members().size();}}

  // Expansion
  clock_t start = clock();
  int curr_vertex = 0;
  for(std::vector< Node >::iterator root_it = root_.begin();
      root_it != root_.end(); ++root_it, ++curr_vertex)
    {
      if(root_it->has_children(curr_vertex)) 
  { siblings_expansion(root_it->children(), dim_max-1); }
    }

  clock_t end = clock();
  std::cout << "Computational time for Rips construction = " << 
    (double)(end - start)/(double)CLOCKS_PER_SEC << std::endl;
}

template<class N>
void 
Flag_simplex_tree<N>::
siblings_expansion(Siblings * siblings, //must contain elements
       int k)
{
  //if (k==0 || members_.empty()) return;
  if(k == 0) return;

  Dictionary_it next = siblings->members().begin(); ++next;

  static std::vector< std::pair<Vertex , Node> > inter;

  for(Dictionary_it s_h = siblings->members().begin();
      s_h != siblings->members().end(); ++s_h,++next)
    {
      if(root_[s_h->first].has_children(s_h->first))
  {
    intersection(inter,  //output intersection
           next,//begin
           siblings->members().end(),//end
           root_[s_h->first].children()->members().begin(),
           root_[s_h->first].children()->members().end(),
           s_h->second.filtration());
    if(inter.size() != 0)
      {
        size_cpx_ += inter.size();
        Siblings * new_sib = new Siblings(siblings,//oncles
            s_h->first, //parent
            inter);//boost::container::ordered_unique_range_t
        inter.clear();
        s_h->second.assign_children(new_sib);
        siblings_expansion(new_sib,k-1);
      }
    else {
      s_h->second.assign_children(siblings); //ensure the children property
      inter.clear();
    }
  }
    }
}

template<class N>
void 
Flag_simplex_tree<N>::
intersection(std::vector< std::pair< typename Flag_simplex_tree<N>::Vertex,
                     typename Flag_simplex_tree<N>::Node > > &       intersection,
             typename Flag_simplex_tree<N>::Dictionary_it            begin1,
             typename Flag_simplex_tree<N>::Dictionary_it            end1,
             typename Flag_simplex_tree<N>::Dictionary_it            begin2,
             typename Flag_simplex_tree<N>::Dictionary_it            end2,
             typename Flag_simplex_tree<N>::Filtration_value         filtration)
{
  if(begin1 == end1 || begin2 == end2) return;// 0;
  while( true )
    {
      if( begin1->first == begin2->first )
  {
    intersection.push_back(std::pair< Vertex, Node >(begin1->first,
                 Node(maximum(begin1->second.filtration(),
                  begin2->second.filtration(),
                  filtration))));
    ++begin1;
    ++begin2;
    if( begin1 == end1 || begin2 == end2 ) return;
  }
      else { 
  if( begin1->first < begin2->first ) 
    {
      ++begin1;
      if(begin1 == end1) return;
    }
  else {
    ++begin2;
    if(begin2 == end2) return;
  }
      }
    }
}

template<class N>
bool 
Flag_simplex_tree<N>::
compare_simplices_fil(const Simplex_handle sh1,
          const Simplex_handle sh2)
{
  if(sh1->second.filtration() != sh2->second.filtration())
    {return sh1->second.filtration() < sh2->second.filtration();}

  return !(is_subface(sh1,sh2)); //is sh1 a subface of sh2
}

template<class N>
bool 
Flag_simplex_tree<N>::
is_subface(Simplex_handle sh1, Simplex_handle sh2)
{
  Simplex_vertex_range rg1 = simplex_vertex_range(sh1);
  Simplex_vertex_range rg2 = simplex_vertex_range(sh2);

  Simplex_vertex_iterator it1 = rg1.begin();
  Simplex_vertex_iterator it2 = rg2.begin();
  while(it1 != rg1.end() && it2 != rg2.end())
    {
      if(*it1 < *it2) {++it2;}
      else {
  if(*it1 == *it2) {++it1; ++it2;}
  else {return false;}
      }
    }
  if(it1 == rg1.end()) return true;
  return false;
}

template<class N>
void 
Flag_simplex_tree<N>::
initialize_filtration()
{
  filtration_vect_.reserve(size_cpx_ - nb_vertices_);      //Attention gestion vertices
  Complex_simplex_range rg = Complex_simplex_range();
  for(Complex_simplex_iterator it = rg.begin();
      it != rg.end(); ++it)
    { filtration_vect_.push_back(*it);}
  stable_sort(filtration_vect_.begin(),filtration_vect_.end(),compare_simplices_fil);
}




