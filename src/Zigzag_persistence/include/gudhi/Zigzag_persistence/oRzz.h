
#include <iostream>
#include <boost/progress.hpp>

#include "io.h"
#include "Simplex_tree.h"

void furthest_point_sampling( std::vector< std::vector<double> > &points
                            , std::vector< std::vector<double> > &f_points
                            , std::vector< double > &epsilon_val);

double dist_sq(std::vector<double> &p, std::vector<double> &q)
{
  double dist=0.; int dim = p.size(); double temp;
  for(int i=0;i<dim;i++) { temp = p[i]-q[i]; dist += temp*temp; }
  return dist;
}


void unique(std::vector < std::vector<double> > & p)
{
  auto cmp_vect = [](std::vector<double> u,std::vector<double> v)->bool
                    {for(int i = 0; i<u.size(); ++i) 
                         { if (u[i]!=v[i]) return u[i]!=v[i]; }
                     return false; };

  std::set < std::vector<double>, decltype(cmp_vect) > s(cmp_vect);
  for(auto &v : p) s.insert(v);
  std::cout << "Unique: Size p = " << p.size() << "   Size s = " << s.size() <<std::endl;
}


// void proc_prob_slow( double threshold
//                    , std::vector< Simplex_tree<>::Simplex_handle > & simpl
//                    , Simplex_tree<> & st
//                    , int &tmp_num) {
//       for(auto sh : st.complex_simplex_range()) 
//       {
//         ++tmp_num;
//         if(st.filtration(sh) > threshold) { simpl.push_back(sh); }
//       }
// }

// void subroutine_prob_slow( std::list< Simplex_tree<>::Simplex_handle > * h_l
//                          , std::vector<Simplex_tree<>::Simplex_handle>::reverse_iterator rit
//                          , Simplex_tree<> &st
//                          , std::list< Simplex_tree<>::Simplex_handle >::iterator & l_itbis)
// {
//   auto l_it = h_l->begin();
//   for(; st.key(*l_it) != st.key(*rit); ++l_it) {}
//   l_itbis = l_it;
// }

// void func_slow_maybe( std::set < Simplex_tree<>::Simplex_handle, bool (*)(Simplex_tree<>::Simplex_handle, Simplex_tree<>::Simplex_handle) >::iterator & l_it
//                     , std::set < Simplex_tree<>::Simplex_handle, bool (*)(Simplex_tree<>::Simplex_handle, Simplex_tree<>::Simplex_handle) >* h_l
//                     , std::vector< Simplex_tree<>::Simplex_handle >::reverse_iterator &rit
//                     , Simplex_tree<> & st )
// {

//    l_it = h_l->begin();
// //   for(; st.key(*l_it) != st.key(*rit); ++l_it) {}
// //   for(int i=0;i<10000000;i++){std::cout << i;}
//   l_it = h_l->find(*rit);
//    return;
// }

// void proc_prob_slow_two(std::vector< Simplex_tree<>::Simplex_handle > & simpl
//                        , Simplex_tree<> & st
//                        , std::ofstream & out
//                        , int & next_key)
// {
//   for(auto rit = simpl.rbegin(); rit != simpl.rend(); ) 
//   {
// //        out << "-1 " << st.key(*rit) << " " << sqrt(rho_sq * epsilon_val[i]) << " \n";
// //    out << "-1 " << st.key(*rit) << " " << sqrt(rho_sq * epsilon_val[i-1]) << " \n";
// //        out << "-1 " << st.key(*rit) << " " << i << " \n";
//     ++next_key; --(st.num_simplices_);

//     //remove from h_list
//     auto h_l = st.node_h_links_[(*rit)->first];
 


//     auto l_it = h_l->begin();
//   l_it = h_l->find(*rit);
//   //  for(; st.key(*l_it) != st.key(*rit); ++l_it) {}
 



//     h_l->erase(l_it);

//     auto tmp_sh = *rit;         ++rit;
//     auto curr_sib = st.self_siblings(tmp_sh);
//     curr_sib->members_.erase(tmp_sh);

//     if(curr_sib->members_.empty() && curr_sib->oncles_ != NULL) 
//     {
//       curr_sib->oncles_->members_[curr_sib->parent_].children_ = curr_sib->oncles_;
//       delete curr_sib;
//     }
//   }
// }



void oRzz(std::string pointfile, 
          std::string outfile,
          double eta_,
          double rho_,
          int dim_max_) 
{
//  std::string pointfile = "../examples/Kl.txt"; //no doublon in file

//  std::string pointfile = "../examples/spiral_4d_2k.dxyz";//../examples/Kl.txt"; //no doublon in file
// std::string pointfile = "../examples/dataset_r50k_15nn_30pc_r600_unique.xyz";

  int dim_max = dim_max_;
  double eta = eta_;   //eta <= rho
  double rho = rho_;

//  std::string outfile = "../examples/Hasse_spiral.txt";
  std::ofstream out(outfile,std::ios::trunc);

  std::vector< std::vector<double> > points;
  read_points(pointfile,points);

  unique(points);

  std::vector<double>                epsilon_val; epsilon_val.reserve(points.size());
  std::vector< std::vector<double> > f_points;    f_points.reserve(points.size());
  
  furthest_point_sampling(points, f_points, epsilon_val);

/////////
//   int i=0;
//   for(auto p : f_points) {
//    for(auto c_p : p) { std::cout << c_p << " "; }
//      std::cout << std::endl; ++i;
// //     std::cout << "     " << sqrt(epsilon_val[i]) << std::endl; ++i;
//   }
//   std::cout << std::endl;  //   epsilon_val[f_points.size()-1] = 0.048; 
////////


  double eta_sq = eta*eta;
  double rho_sq = rho*rho;
  int next_key  = 0;

  Simplex_tree <> st;
  //insert all vertices
  for(int i = 0; i < f_points.size(); ++i) 
  { 
    out << "0 " << next_key << " " << sqrt(epsilon_val[i]) << " \n"; //filtration func based
//    out << "0 " << next_key << " " << i << " \n"; //index based...

    auto sh = st.oRzz_insert_vertex(i , 0); //fil = 0 to not get remove in the oRzz
    st.assign_key(sh,next_key);
    ++next_key;
  } //st.oRzz_insert_vertex (0, next_key, out);


// NB: index zz_persistence:    A -> B <- C -> D <- E -> F ...
//                              0  [1;1]     [2;2]     [3;3]
// i.e everything added in A->B has filtration i
//     everything removed in B <- C has filtration i  

  (st.num_simplices_) = 1; //one vertex
  int max_num_simplices = 1;//count only one vertex (0) for start  //f_points.size();
  std::vector< Simplex_tree<>::Simplex_handle > simpl;
  for(int i = 1; i < f_points.size(); ++i) 
  {
    //Up right arrow
    ++(st.num_simplices_); //new vertex i
    //all edges with i
    for(int v = 0; v < i; ++v) 
    {
      double dist_uv = dist_sq(f_points[i],f_points[v]);
      if( dist_uv <= rho_sq * epsilon_val[i-1] )//&&(dist_uv >  eta_sq * epsilon_val[i-1])) 
      { //new edge //not epsilon i ?
        st.oRzz_insert_edge_expand( v, i, dist_uv, dim_max, simpl );    //in tree distances are squared
      }
    }

    //between already existing vertices, with the bigger radius
    for(int v = 0; v < i; ++v) 
    {
      for(int u = 0; u < v; ++u) 
      {
        double dist_uv = dist_sq(f_points[u],f_points[v]);
        if( (dist_uv <= rho_sq * epsilon_val[i-1]) && (dist_uv > eta_sq * epsilon_val[i-1])) 
        {  
          st.oRzz_insert_edge_expand(u, v, dist_uv, dim_max, simpl); 
        }
      }
    }

    //sort and update the keys now + write in file
//    auto beg = st.filtration_vect_.begin()+num_simp_before_insert;
  
    // std::list<Simplex_tree<>::Simplex_handle> tmp_tail;
    // tmp_tail.splice(tmp_tail.end(), st.filtration_vect_, beg, st.filtration_vect_.end() ); //cut the tail
    // tmp_tail.sort(Simplex_tree<>::is_before_in_filtration(&st));

    stable_sort(simpl.begin(),simpl.end(),Simplex_tree<>::is_before_in_filtration_zz(&st));
    
    // auto tmp_sh_it = simpl.begin();
    // for(auto sh_it = simpl.begin(); ;) {
    //   tmp_sh_it = sh_it;
    //   ++sh_it;
    //   if(sh_it == simpl.end()) { break; }
    //   if( Simplex_tree<>::is_before_in_filtration_zz(&st)(tmp_sh_it, sh_it) ) 
    // }



    for(auto sh : simpl)
    { 
      st.assign_key(sh, next_key); 

      out << st.dimension(sh) << " ";
      for(auto b_sh : st.boundary_simplex_range(sh)) { out << st.key(b_sh) << " "; }

      out << st.key(sh) << " " << sqrt(rho_sq * epsilon_val[i-1]) << " \n"; //filtration based
//      out << st.key(sh) << " " << i << " \n"; //filtration value = index i
      
      //out << st.key(sh) << " " << sqrt(rho_sq * epsilon_val[i-1]) << " \n"; 
      //log2(rho_sq * epsilon_val[i-1])/2. << " \n";//      " << st.filtration(*it) << "\n";
      ++next_key;
    }
//    std::stable_sort(beg,st.filtration_vect_.end(),Simplex_tree<>::is_before_in_filtration(&st));
    // for(auto it = tmp_tail.begin(); it != tmp_tail.end(); ++it)
    // { 
    //   st.assign_key(*it, next_key); 

    //   out << st.dimension(*it) << " ";
    //   for(auto b_sh : st.boundary_simplex_range(*it)) { out << st.key(b_sh) << " "; }
    //   out << st.key(*it) << "       " << st.filtration(*it) << "\n";

    //   ++next_key;
    // }
    //reattach tail
  //  st.filtration_vect_.splice(st.filtration_vect_.end(),tmp_tail);



    std::cout << i << "   " << st.num_simplices() << "    ---    ";

    if(st.num_simplices() > max_num_simplices) { max_num_simplices = st.num_simplices(); }

    // std::cout << "------------------After all inserts: beg \n";
    // std::cout << st ;
    // for(auto sh : st.filtration_simplex_range()) {st.display_simplex(sh); std::cout << std::endl;}
    // std::cout << "------------------end \n";


    //remove everything under threshold
//    std::cout << " +++ rm with eps = " << eta_sq * epsilon_val[i] << " \n";
    simpl.clear();
    double threshold = eta_sq * epsilon_val[i];

    int tmp_num = 0;

    if(threshold < rho_sq * epsilon_val[i-1]) 
    {


//      proc_prob_slow(threshold,simpl,st,tmp_num);
      for(auto sh : st.complex_simplex_range()) 
      {
        ++tmp_num;
        if(st.filtration(sh) > threshold) { simpl.push_back(sh); }
      }




//      if(tmp_num < st.num_simplices()) {std::cout << "Forget simplices \n";}
    


      sort(simpl.begin(),simpl.end(),Simplex_tree<>::is_before_in_filtration_zz(&st));
   



//       proc_prob_slow_two(simpl,st,out,next_key);
      for(auto rit = simpl.rbegin(); rit != simpl.rend(); ) 
      {
//        out << "-1 " << st.key(*rit) << " " << sqrt(rho_sq * epsilon_val[i]) << " \n";
        out << "-1 " << st.key(*rit) << " " << sqrt(rho_sq * epsilon_val[i-1]) << " \n";
//        out << "-1 " << st.key(*rit) << " " << i << " \n";
        ++next_key; --(st.num_simplices_);

        //remove from h_list
        auto * h_l = st.node_h_links_[(*rit)->first];


        // auto l_it = h_l->begin();
        // subroutine_prob_slow(h_l,rit,st,l_it);


//        auto l_it = h_l->begin();
//        for(; st.key(*l_it) != st.key(*rit); ++l_it) {}



        auto l_it = h_l->begin();
        l_it = h_l->find(*rit);

      //  func_slow_maybe(l_it,h_l,rit,st);
      //        auto l_it = h_l->find(*rit);

        h_l->erase(l_it);

        auto tmp_sh = *rit;         ++rit;
        auto curr_sib = st.self_siblings(tmp_sh);
        curr_sib->members_.erase(tmp_sh);

        if(curr_sib->members_.empty() && curr_sib->oncles_ != NULL) 
        {
          curr_sib->oncles_->members_[curr_sib->parent_].children_ = curr_sib->oncles_;
          delete curr_sib;
        }
      }




      simpl.clear();
    



    }

      //st.remove_over_threshold(eta_sq * epsilon_val[i], out, simpl);

      std::cout << st.num_simplices() << "              " << next_key << " \n";
    
  //    std::cout << "------------------After all removes: beg \n";
  //    std::cout << st ;
  //    for(auto sh : st.filtration_simplex_range()) {st.display_simplex(sh); std::cout << std::endl;}
  //    std::cout << "------------------end \n";
  }

  //remove all vertices
  for(int i = f_points.size()-1; i > -1 ; --i) 
  {
    out << "-1 " << i << " " << sqrt(epsilon_val[f_points.size()-1]) << " \n"; //index based...
//    out << "-1 " << i << " " << f_points.size()-1 << " \n";
    ++next_key;
  } 
  out.close();
  std::cout << "Max size complex = " << max_num_simplices << "     next_key = " << next_key << std::endl;
}


void furthest_point_sampling( std::vector< std::vector<double> > &points
                            , std::vector< std::vector<double> > &f_points
                            , std::vector< double > &epsilon_val)
{
  std::vector< double > dist_to_P(points.size(),INFINITY);

  int idx_first_point = 0;
  f_points.push_back(points[idx_first_point]);

  double max_dist;  double curr_dist; int next_point;

  boost::progress_display show_progress(points.size()-1);
  for(int i = 1 ; i < points.size(); ++i) //pick the ith element of f_points
  {
    max_dist = -1;
    
    for(int j = 0; j < points.size(); ++j) 
    {
      curr_dist = dist_sq(f_points[i-1], points[j]);
      if(dist_to_P[j] > curr_dist) { dist_to_P[j] = curr_dist; }
      if(dist_to_P[j] > max_dist)  { max_dist = dist_to_P[j]; next_point = j;}
    }
  
    epsilon_val.push_back( max_dist ); //i-1
    f_points.push_back(points[next_point]);//i
  
    ++show_progress;
  }
  epsilon_val.push_back(0);
}
