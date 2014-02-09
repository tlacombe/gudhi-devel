/*
 *  Simplex_tree.hpp
 *  Gudhi
 *
 *  Created by Clément Maria on 1/7/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */


// Print a Siblings in os, recursively.
template <class V, class F, class N, class MC > 
std::ostream& operator<<(std::ostream& os, 
                         Simplex_tree_siblings<V,F,N,MC> & obj)
{
  os << "--Oncles: @ " << (long int)(obj.oncles()) << "\n";
  os << "--Parent:   " << obj.parent() << "\n";
  os << "Siblings: @ " << (long int)(&obj) << "\n";
  for(typename MC::iterator sh = obj.members().begin();
      sh != obj.members().end(); ++sh)
    {    os << "[" << sh->first << ":" << sh->second.filtration() <<"] ";    }
  os << std::endl << std::endl;
  for(typename MC::iterator sh = obj.members().begin();
      sh != obj.members().end(); ++sh)
    {if(sh->second.has_children(sh->first)) os << *(sh->second.children());}
  return os;
}
// Print a Simplex_tree in os.
template< class MS >
std::ostream& operator<<(std::ostream& os, 
                         Simplex_tree< MS > & obj)
{
  os << "Simplex Tree: \n";
  os << "Size Cpx   = " << obj.size_complex() << std::endl;
  os << "nb_V       = " << obj.nb_vertices() << std::endl;
  os << std::endl;
  os << obj.root();
return os;
}    













