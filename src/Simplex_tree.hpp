/*
 *  Simplex_tree.hpp
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

template<class MS>
bool 
Simplex_tree<MS>::
is_subface(Simplex_handle sh1, Simplex_handle sh2)
{
  Simplex_vertex_range< Simplex_tree >    rg1 = simplex_vertex_range(sh1);
  Simplex_vertex_range< Simplex_tree >    rg2 = simplex_vertex_range(sh2);
  Simplex_vertex_iterator< Simplex_tree > it1 = rg1.begin();
  Simplex_vertex_iterator< Simplex_tree > it2 = rg2.begin();
  while(it1 != rg1.end() && it2 != rg2.end()) {
      if(*it1 < *it2) {++it2;}
      else {
          if(*it1 == *it2) {++it1; ++it2;}
          else {return false;}
      }
    }
  if(it1 == rg1.end()) return true;
  return false;
}

/**
* \brief Inserts a Simplex_handle for every simplex in the simplicial
* complex, except vertices, and sort they according to the filtration.
*
* The use of a depth-first traversal of the simplex tree combined with
* a stable sort is meant to optimize the order of simplices with same
* filtration value. The heuristic consists in inserting the cofaces of a
* simplex as soon as possible.
*
* \todo Check if the heuristic is efficient. How should we deal with vertices? 
*
*/
template<class MS>
void 
Simplex_tree<MS>::
initialize_filtration()
{
  filtration_vect_.reserve(size_cpx_ - nb_vertices_);      //Attention gestion vertices
  Complex_simplex_range< Simplex_tree > rg = Complex_simplex_range< Simplex_tree >();
  for(Complex_simplex_iterator< Simplex_tree > it = rg.begin();
      it != rg.end(); ++it)
    { filtration_vect_.push_back(*it);}
  stable_sort(filtration_vect_.begin(),filtration_vect_.end(),is_before_in_filtration(this));
}



/**
* Use of "static" makes the librarynot thread safe.
*
* \todo non-static?
*/
template < class MS >
void 
Simplex_tree<MS>::
siblings_expansion(Siblings * siblings, //must contain elements
                   int k)
{
  //if (k==0 || members_.empty()) return;
  if(k == 0) return;
  Dictionary_it next = siblings->members().begin(); ++next;

  static std::vector< std::pair<Vertex , Node> > inter; // <-------static

  for(Dictionary_it s_h = siblings->members().begin();
      s_h != siblings->members().end(); ++s_h,++next)
    {
      if(this->root_[s_h->first].has_children(s_h->first))
      {
        intersection(inter,  //output intersection
                     next,                     //begin
                     siblings->members().end(),//end
                     this->root_[s_h->first].children()->members().begin(),
                     this->root_[s_h->first].children()->members().end(),
                     s_h->second.filtration());
        if(inter.size() != 0)
        {
          this->size_cpx_ += inter.size();
          Siblings * new_sib = new Siblings(siblings,   //oncles
                                            s_h->first, //parent
                                            inter);     //boost::container::ordered_unique_range_t
          inter.clear();
          s_h->second.assign_children(new_sib);
          siblings_expansion(new_sib,k-1);
        }
        else 
        { s_h->second.assign_children(siblings); //ensure the children property
          inter.clear();}
      }
    }
}


/**
* \brief Print a Siblings in os.
*/
template <class V, class F, class N> 
std::ostream& operator<<(std::ostream& os, 
                         Simplex_tree_siblings<V,F,N> & obj)
{
  os << "--Oncles: @ " << (long int)(obj.oncles()) << "\n";
  os << "--Parent:   " << obj.parent() << "\n";
  os << "Siblings: @ " << (long int)(&obj) << "\n";
  for(typename Simplex_tree_siblings<V,F,N>::Dictionary::iterator sh = obj.members().begin();
      sh != obj.members().end(); ++sh)
    {    os << "[" << sh->first << ":" << sh->second.filtration() <<"] ";    }
  os << std::endl << std::endl;
  for(typename Simplex_tree_siblings<V,F,N>::Dictionary::iterator sh = obj.members().begin();
      sh != obj.members().end(); ++sh)
    {if(sh->second.has_children(sh->first)) os << *(sh->second.children());}
  return os;
}
/**
* Print a Simplex_tree in os.
*/
template< class MS >
std::ostream& operator<<(std::ostream& os, 
                         Simplex_tree< MS > & obj)
{
  os << "Simplex Tree: \n";
  os << "Size Cpx   = " << obj.size_complex() << std::endl;
  os << "nb_V    = " << obj.nb_vertices() << std::endl;
  os << std::endl;

  int v = 0;
  os << "@ ROOT:   ";
  for(typename std::vector< typename Simplex_tree< MS >::Node >::iterator it = obj.root().begin();
      it != obj.root().end(); ++it, ++v)
    {    os << v << " ";  }
  os << std::endl << std::endl;;

  v = 0;
  for(typename std::vector< typename Simplex_tree< MS >::Node >::iterator it = obj.root().begin();
      it != obj.root().end(); ++it,++v)
    {
      if(it->has_children(v)) os << *(it->children());    
    }
  return os;
}    













