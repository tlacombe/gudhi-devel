
#include <iostream>
#include <fstream>

#ifdef GUDHI_USE_TBB
#include <tbb/tbb.h>
#endif

/** Given a set of points p_1, ... , p_n ordered by their insertion order in the 
  * filtration, 
  * computes the edge-filtration corresponding to the oscillating Rips zigzag 
  * filtration of the set of points, i.e.,
  * ... <- R({p_0, ... , p_{i}}, nu * eps_i) ->
  *                   R({p_0, ... , p_i, p_{i+1}}, mu * eps_i) <- 
  *                              R({p_0, ... , p_i, p_{i+1}}, nu * eps_{i+1}) -> ...
  * where 0 < nu <= mu, and eps_i is defined as the sparsity of the point cloud 
  * {p_1, ... , p_i}, i.e., the shortest distance between two points in the set. 
  * This is a decreasing sequence of numbers.
  *
  * The function computes the eps_i in filtration_value[i] = eps_i, with eps_0 = 
  * infinity. A simplex appearing in the inclusion  
  * R({p_0, ... , p_{i}}, nu * eps_i) -> R({p_1, ... , p_i, p_{i+1}}, mu * eps_i)
  * is given filtration value eps_i.
  *
  * filtration_values must be empty, receives the filtration values.
  * edge_filtration must be empty, receives the edge filtration.
  */
template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t >
void points_to_edge_filtration(Point_container const        &points,
                               Distance                      distance,
                               double                        nu,
                               double                        mu,
                               std::vector<FiltrationValue> &filtration_values,
                               std::vector<Edge_t>          &edge_filtration )
{
  //computes the eps_i = min_{p,q \in P_i} d(p,q) naively, in parallel
  size_t n = points.size();
  filtration_values.resize(n);
  filtration_values[0] = std::numeric_limits<double>::infinity();//eps_0
  for(size_t i = 1; i < n; ++i) {//truly parallelisable
    double dist = std::numeric_limits<double>::infinity();
    for(size_t j = 0; j < i; ++j) { //d(P_{i-1}, p_i)
      auto curr_dist = distance(points[i],points[j]); 
      if(dist > curr_dist) { dist = curr_dist; } //maintain shortest distance
    }
    filtration_values[i] = dist; //d(P_{i-1}, p_i)
  }
  //turn filtration_value[i] into sparsity of {p_0, ... , p_i}
  for(size_t i = 1; i < n; ++i) {
    if(filtration_values[i] > filtration_values[i-1]) //make decreasing
    {  filtration_values[i] = filtration_values[i-1];  }
  }
  //initialise R({p_0}, \nu * eps_0)
  edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add p_0
  
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_1, ... , p_i}, mu * eps_i)   radius   
    for(size_t j = 1; j <= i; ++j) {//nu eps_i < length(p_j, p_k) <= mu eps_i
      for(size_t k = 0; k < j; ++k) {
        if(distance(points[j],points[k]) <= mu * filtration_values[i] && 
           distance(points[j],points[k]) > nu * filtration_values[i]) {
          edge_filtration.emplace_back(k, j, filtration_values[i], true);//edge k,j 
        }
      }
    }

    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_1, ... , p_i, p_i+1}, mu * eps_i)
    edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(size_t j = 0; j < i+1; ++j) {//edges (p_{i+1}, p_j) of length <= mu * eps_i 
      if(distance(points[j],points[i+1]) <= mu * filtration_values[i]) {
        edge_filtration.emplace_back(j, i+1, filtration_values[i], true);//edge 
      }
    }
    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    // std::cout << "Remove edges going down eps_" << i << "to eps_" << i+1 <<" \n";
    for(size_t j = 1; j <= i+1; ++j) {//nu eps_i+1 < length(p_j, p_k) <= mu eps_i
      for(size_t k = 0; k < j; ++k) {
        auto dist = distance(points[j],points[k]);
        if(dist <= mu *filtration_values[i] && dist > nu *filtration_values[i+1]) {
          // std::cout << "  " << j << " " << k << "   " << nu * filtration_values[i+1] << " < *" << dist << "* <= " << mu * filtration_values[i] << "\n";
          edge_filtration.emplace_back(k, j, filtration_values[i+1], false);
        }
      }
    }
    // std::cout << "done.\n";
  }
}

template<typename FiltrationValue>
struct point_distance_cmp {
  bool operator()( std::pair<int, FiltrationValue> p1
                 , std::pair<int, FiltrationValue> p2 ) {
    { 
      if(p1.second != p2.second) {return p1.second < p2.second;} //shorter first
      return p1.first < p2.first; 
    }
  }
};

template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t > 
struct edge_cmp {

  edge_cmp(Point_container &points, Distance distance) 
  : points_(&points), distance_(distance) {}

  bool operator()(Edge_t e1, Edge_t e2) 
  {
    FiltrationValue dist1 = distance_((*points_)[e1.u()], (*points_)[e1.v()]);//l(e1)  
    FiltrationValue dist2 = distance_((*points_)[e2.u()], (*points_)[e2.v()]);//l(e2)
    if(dist1  != dist2)  {return dist1 < dist2;}
    if(e1.u() != e2.u()) {return e1.u() < e2.u();}
    return e1.v() < e2.v(); 
  }

private:
  Point_container * points_;
  Distance distance_;
};


#ifdef GUDHI_USE_TBB //parallel version
template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t >
void real_fast_points_to_edge_filtration(Point_container         &points,
                               Distance                      distance,
                               double                        nu,
                               double                        mu,
                               std::vector<FiltrationValue> &filtration_values,
                               std::vector<Edge_t>          &edge_filtration )
{
  //computes the eps_i = min_{p,q \in P_i} d(p,q) naively, in parallel
  size_t n = points.size();
  filtration_values.resize(n);
  filtration_values[0] = std::numeric_limits<double>::infinity();//eps_0
  point_distance_cmp<FiltrationValue> cmp;
  size_t number_of_arrows = 0;

  //syntax tbb
  // for( size_t i=0; i!=n; ++i ) {Foo(a[i]);}
  // tbb::parallel_for( size_t(0), n, [&]( size_t i ) {Foo(a[i]);} );

  // for(int i = 1; i < n; ++i) {//truly parallelisable
  tbb::parallel_for(size_t(1), n, [&](size_t i) {//O(n^2)
    double dist = std::numeric_limits<double>::infinity();
    for(size_t j = 0; j < i; ++j) { //d(P_{i-1}, p_i)
      auto curr_dist = distance(points[i],points[j]); 
      if(dist > curr_dist) { dist = curr_dist; } //maintain shortest distance
    }
    filtration_values[i] = dist; //d(P_{i-1}, p_i)
  });

  //turn filtration_value[i] into sparsity of {p_0, ... , p_i}
  for(size_t i = 1; i < n; ++i) {
    if(filtration_values[i] > filtration_values[i-1]) //make decreasing
    {  filtration_values[i] = filtration_values[i-1];  }
  }

  //initialize distance matrix
  //dist_matrix[i] contains all pairs (j, d(p_i,p_j)) for j < i
  std::vector< std::vector< std::pair<int, FiltrationValue> > > dist_matrix;
  dist_matrix.resize(n);
  for(size_t i = 0; i < n; ++i) {//for all vertices
    dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
  }

  // for(size_t i = 0; i < n; ++i) {//all vertices
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    for(size_t j = 0; j < i; ++j) {
      dist_matrix[i][j] 
                = std::pair<int, FiltrationValue>(j, distance(points[j],points[i]));
    }
    std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp);
  } );

  typename std::vector< std::pair<int, FiltrationValue> >::iterator it;

//edges[i]==list of edges added and removed at eps_i
  std::vector< std::vector< Edge_t > > edges_added; 
  std::vector< std::vector< Edge_t > > edges_removed;
  edges_added.resize(n);  edges_removed.resize(n);

  // size_t m = n*(n-1);
  // for(size_t i=0; i<n; ++i) {
  //   edges_added[i].reserve(m);     edges_removed[i].reserve(m);
  // }

  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
  // tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_1, ... , p_i}, mu * eps_i)   radius   
    for(size_t j = 1; j <= i ; ++j) {
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
                std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp);//first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i]) { break; }//short edge there
        // edge_filtration.emplace_back(it->first, j, filtration_values[i], true);
        edges_added[i].emplace_back(it->first, j, filtration_values[i], true);
        ++number_of_arrows;
      }
    }
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_1, ... , p_i, p_i+1}, mu * eps_i)
    // edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp); //first striclty longer edge
    while(it != dist_matrix[i+1].begin()) {
      --it;
      // edge_filtration.emplace_back(it->first, i+1, filtration_values[i], true);
      edges_added[i].emplace_back(it->first, i+1, filtration_values[i], true);
      ++number_of_arrows;
    }

    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    for(size_t j = 1; j <= i+1; ++j) {//nu eps_i+1 < length(p_j, p_k) <= mu eps_i
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp ); //first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i+1]) { break; }
        // edge_filtration.emplace_back(it->first, j, filtration_values[i+1], false);
        edges_removed[i].emplace_back(it->first, j, filtration_values[i+1], false);
        ++number_of_arrows;
      }
    }
  } 
  // );

//sort insertions and deletion by edge length
  edge_cmp<Point_container, Distance, FiltrationValue, Edge_t > e_cmp(points, distance);
  // for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
  tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
    //add shortest edges first
    std::stable_sort(edges_added[i].begin(), edges_added[i].end(), e_cmp);
    //remove longest edges first (read from right to left)
    std::stable_sort(edges_removed[i].begin(), edges_removed[i].end(), e_cmp);
  } );

  number_of_arrows += n; //count vertices additions
  edge_filtration.reserve(number_of_arrows);

  //initialise R({p_0}, \nu * eps_0)
  edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add p_0
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(auto edg_it = edges_added[i].begin(); 
            edg_it != edges_added[i].end(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
    for(auto edg_it = edges_removed[i].rbegin(); 
             edg_it != edges_removed[i].rend(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
  }
  std::cout << "done.\n";
}
#else //non parallel version
template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t >
void real_fast_points_to_edge_filtration(Point_container         &points,
                               Distance                      distance,
                               double                        nu,
                               double                        mu,
                               std::vector<FiltrationValue> &filtration_values,
                               std::vector<Edge_t>          &edge_filtration )
{
  //computes the eps_i = min_{p,q \in P_i} d(p,q) naively, in parallel
  size_t n = points.size();
  filtration_values.resize(n);
  filtration_values[0] = std::numeric_limits<double>::infinity();//eps_0
  point_distance_cmp<FiltrationValue> cmp;
  size_t number_of_arrows = 0;

  for(size_t i = 1; i < n; ++i) {
    double dist = std::numeric_limits<double>::infinity();
    for(size_t j = 0; j < i; ++j) { //d(P_{i-1}, p_i)
      auto curr_dist = distance(points[i],points[j]); 
      if(dist > curr_dist) { dist = curr_dist; } //maintain shortest distance
    }
    filtration_values[i] = dist; //d(P_{i-1}, p_i)
  }
  //turn filtration_value[i] into sparsity of {p_0, ... , p_i}
  for(size_t i = 1; i < n; ++i) {
    if(filtration_values[i] > filtration_values[i-1]) //make decreasing
    {  filtration_values[i] = filtration_values[i-1];  }
  }
  //dist_matrix[i] contains all pairs (j, d(p_i,p_j)) for j < i
  std::vector< std::vector< std::pair<int, FiltrationValue> > > dist_matrix;
  dist_matrix.resize(n);
  for(size_t i = 0; i < n; ++i) {//for all vertices
    dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
  }
  //sort neighbours by increasing distance
  for(size_t i = 0; i < n; ++i) {//all vertices
    for(size_t j = 0; j < i; ++j) {
      dist_matrix[i][j] = 
                  std::pair<int, FiltrationValue>(j, distance(points[j],points[i]));
    }
  std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp );
  }

  typename std::vector< std::pair<int, FiltrationValue> >::iterator it;

//edges[i]==list of edges added and removed at eps_i
  std::vector< std::vector< Edge_t > > edges_added; 
  std::vector< std::vector< Edge_t > > edges_removed;
  edges_added.resize(n);  edges_removed.resize(n);

  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_1, ... , p_i}, mu * eps_i)   radius   
    for(size_t j = 1; j <= i ; ++j) {
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
                std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp);//first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i]) { break; }//short edge there
        // edge_filtration.emplace_back(it->first, j, filtration_values[i], true);
        edges_added[i].emplace_back(it->first, j, filtration_values[i], true);
        ++number_of_arrows;
      }
    }
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_1, ... , p_i, p_i+1}, mu * eps_i)
    // edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp); //first striclty longer edge
    while(it != dist_matrix[i+1].begin()) {
      --it;
      // edge_filtration.emplace_back(it->first, i+1, filtration_values[i], true);
      edges_added[i].emplace_back(it->first, i+1, filtration_values[i], true);
      ++number_of_arrows;
    }

    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    for(size_t j = 1; j <= i+1; ++j) {//nu eps_i+1 < length(p_j, p_k) <= mu eps_i
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp ); //first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i+1]) { break; }
        // edge_filtration.emplace_back(it->first, j, filtration_values[i+1], false);
        edges_removed[i].emplace_back(it->first, j, filtration_values[i+1], false);
        ++number_of_arrows;
      }
    }
  }

  edge_cmp<Point_container, Distance, FiltrationValue, Edge_t > e_cmp(points, distance);

  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    //add shortest edges first
    std::stable_sort(edges_added[i].begin(), edges_added[i].end(), e_cmp);
    //remove longest edges first (read from right to left)
    std::stable_sort(edges_removed[i].begin(), edges_removed[i].end(), e_cmp);
  }

  number_of_arrows += n; //count vertices additions
  edge_filtration.reserve(number_of_arrows);

  //initialise R({p_0}, \nu * eps_0)
  edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add p_0
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(auto edg_it = edges_added[i].begin(); 
            edg_it != edges_added[i].end(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
    for(auto edg_it = edges_removed[i].rbegin(); 
             edg_it != edges_removed[i].rend(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
  }
}
#endif


//puts all vertices first, with distinct filtration value <-- to eventually remove
#ifdef GUDHI_USE_TBB //parallel version
template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t >
void fast_points_to_edge_filtration(Point_container         &points,
                               Distance                      distance,
                               double                        nu,
                               double                        mu,
                               std::vector<FiltrationValue> &filtration_values,
                               std::vector<Edge_t>          &edge_filtration )
{
  //computes the eps_i = min_{p,q \in P_i} d(p,q) naively, in parallel
  size_t n = points.size();
  filtration_values.resize(n);
  filtration_values[0] = std::numeric_limits<double>::max();//std::numeric_limits<double>::infinity();//eps_0
  point_distance_cmp<FiltrationValue> cmp;
  size_t number_of_arrows = 0;

  //syntax tbb
  // for( size_t i=0; i!=n; ++i ) {Foo(a[i]);}
  // tbb::parallel_for( size_t(0), n, [&]( size_t i ) {Foo(a[i]);} );

  // for(int i = 1; i < n; ++i) {//truly parallelisable
  tbb::parallel_for(size_t(1), n, [&](size_t i) {//O(n^2)
    double dist = std::numeric_limits<double>::infinity();
    for(size_t j = 0; j < i; ++j) { //d(P_{i-1}, p_i)
      auto curr_dist = distance(points[i],points[j]); 
      if(dist > curr_dist) { dist = curr_dist; } //maintain shortest distance
    }
    filtration_values[i] = dist; //d(P_{i-1}, p_i)
  });

  //turn filtration_value[i] into sparsity of {p_0, ... , p_i}
  for(size_t i = 1; i < n; ++i) {
    if(filtration_values[i] > filtration_values[i-1]) //make decreasing
    {  filtration_values[i] = filtration_values[i-1];  }
  }

  //initialize distance matrix
  //dist_matrix[i] contains all pairs (j, d(p_i,p_j)) for j < i
  std::vector< std::vector< std::pair<int, FiltrationValue> > > dist_matrix;
  dist_matrix.resize(n);
  for(size_t i = 0; i < n; ++i) {//for all vertices
    dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
  }

  // for(size_t i = 0; i < n; ++i) {//all vertices
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    for(size_t j = 0; j < i; ++j) {
      dist_matrix[i][j] 
                = std::pair<int, FiltrationValue>(j, distance(points[j],points[i]));
    }
    std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp);
  } );

  typename std::vector< std::pair<int, FiltrationValue> >::iterator it;

//edges[i]==list of edges added and removed at eps_i
  std::vector< std::vector< Edge_t > > edges_added; 
  std::vector< std::vector< Edge_t > > edges_removed;
  edges_added.resize(n);  edges_removed.resize(n);

  // size_t m = n*(n-1);
  // for(size_t i=0; i<n; ++i) {
  //   edges_added[i].reserve(m);     edges_removed[i].reserve(m);
  // }

  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
  // tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_0, ... , p_i}, mu * eps_i)   radius   
    for(size_t j = 1; j <= i ; ++j) {
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
                std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp);//first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i]) { break; }//short edge there
        // edge_filtration.emplace_back(it->first, j, filtration_values[i], true);
        edges_added[i].emplace_back(it->first, j, filtration_values[i], true);
        ++number_of_arrows;
      }
    }
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_0, ... , p_i, p_i+1}, mu * eps_i)
    // edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp); //first striclty longer edge
    while(it != dist_matrix[i+1].begin()) {
      --it;
      // edge_filtration.emplace_back(it->first, i+1, filtration_values[i], true);
      edges_added[i].emplace_back(it->first, i+1, filtration_values[i], true);
      ++number_of_arrows;
    }

    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    for(size_t j = 1; j <= i+1; ++j) {//nu eps_i+1 < length(p_j, p_k) <= mu eps_i
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp ); //first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i+1]) { break; }
        // edge_filtration.emplace_back(it->first, j, filtration_values[i+1], false);
        edges_removed[i].emplace_back(it->first, j, filtration_values[i+1], false);
        ++number_of_arrows;
      }
    }
  } 
  // );

//sort insertions and deletion by edge length
  edge_cmp<Point_container, Distance, FiltrationValue, Edge_t > e_cmp(points, distance);
  // for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
  tbb::parallel_for(size_t(0), n-1, [&](size_t i) {
    //add shortest edges first
    std::stable_sort(edges_added[i].begin(), edges_added[i].end(), e_cmp);
    //remove longest edges first (read from right to left)
    std::stable_sort(edges_removed[i].begin(), edges_removed[i].end(), e_cmp);
  } );

  number_of_arrows += n; //count vertices additions
  edge_filtration.reserve(number_of_arrows);

//all vertices
  for(size_t i=0; i<n; ++i) {
    edge_filtration.emplace_back(i, i, std::numeric_limits<double>::infinity(), true);    
  }
  // edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add p_0
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    // edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(auto edg_it = edges_added[i].begin(); 
            edg_it != edges_added[i].end(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
    for(auto edg_it = edges_removed[i].rbegin(); 
             edg_it != edges_removed[i].rend(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
  }
}
#else //non parallel version
template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t >
void fast_points_to_edge_filtration(Point_container         &points,
                               Distance                      distance,
                               double                        nu,
                               double                        mu,
                               std::vector<FiltrationValue> &filtration_values,
                               std::vector<Edge_t>          &edge_filtration )
{
  //computes the eps_i = min_{p,q \in P_i} d(p,q) naively, in parallel
  size_t n = points.size();
  filtration_values.resize(n);
  filtration_values[0] = std::numeric_limits<double>::max();//std::numeric_limits<double>::infinity();//eps_0
  point_distance_cmp<FiltrationValue> cmp;
  size_t number_of_arrows = 0;

  for(size_t i = 1; i < n; ++i) {
    double dist = std::numeric_limits<double>::infinity();
    for(size_t j = 0; j < i; ++j) { //d(P_{i-1}, p_i)
      auto curr_dist = distance(points[i],points[j]); 
      if(dist > curr_dist) { dist = curr_dist; } //maintain shortest distance
    }
    filtration_values[i] = dist; //d(P_{i-1}, p_i)
  }
  //turn filtration_value[i] into sparsity of {p_0, ... , p_i}
  for(size_t i = 1; i < n; ++i) {
    if(filtration_values[i] > filtration_values[i-1]) //make decreasing
    {  filtration_values[i] = filtration_values[i-1];  }
  }
  //dist_matrix[i] contains all pairs (j, d(p_i,p_j)) for j < i
  std::vector< std::vector< std::pair<int, FiltrationValue> > > dist_matrix;
  dist_matrix.resize(n);
  for(size_t i = 0; i < n; ++i) {//for all vertices
    dist_matrix[i] = std::vector< std::pair<int, FiltrationValue> >();
    dist_matrix[i].resize(i);
  }
  //sort neighbours by increasing distance
  for(size_t i = 0; i < n; ++i) {//all vertices
    for(size_t j = 0; j < i; ++j) {
      dist_matrix[i][j] = 
                  std::pair<int, FiltrationValue>(j, distance(points[j],points[i]));
    }
  std::stable_sort(dist_matrix[i].begin(), dist_matrix[i].end(), cmp );
  }

  typename std::vector< std::pair<int, FiltrationValue> >::iterator it;

//edges[i]==list of edges added and removed at eps_i
  std::vector< std::vector< Edge_t > > edges_added; 
  std::vector< std::vector< Edge_t > > edges_removed;
  edges_added.resize(n);  edges_removed.resize(n);

  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_1, ... , p_i}, mu * eps_i)   radius   
    for(size_t j = 1; j <= i ; ++j) {
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
                std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp);//first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i]) { break; }//short edge there
        // edge_filtration.emplace_back(it->first, j, filtration_values[i], true);
        edges_added[i].emplace_back(it->first, j, filtration_values[i], true);
        ++number_of_arrows;
      }
    }
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_1, ... , p_i, p_i+1}, mu * eps_i)
    // edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    it = std::upper_bound(dist_matrix[i+1].begin(), dist_matrix[i+1].end(), 
            std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp); //first striclty longer edge
    while(it != dist_matrix[i+1].begin()) {
      --it;
      // edge_filtration.emplace_back(it->first, i+1, filtration_values[i], true);
      edges_added[i].emplace_back(it->first, i+1, filtration_values[i], true);
      ++number_of_arrows;
    }

    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    for(size_t j = 1; j <= i+1; ++j) {//nu eps_i+1 < length(p_j, p_k) <= mu eps_i
      it = std::upper_bound(dist_matrix[j].begin(), dist_matrix[j].end(), 
             std::pair<int, FiltrationValue>(n, mu * filtration_values[i]), cmp ); //first striclty longer edge
      while(it != dist_matrix[j].begin()) {
        --it;
        if(it->second <= nu * filtration_values[i+1]) { break; }
        // edge_filtration.emplace_back(it->first, j, filtration_values[i+1], false);
        edges_removed[i].emplace_back(it->first, j, filtration_values[i+1], false);
        ++number_of_arrows;
      }
    }
  }

  edge_cmp<Point_container, Distance, FiltrationValue, Edge_t > e_cmp(points, distance);

  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    //add shortest edges first
    std::stable_sort(edges_added[i].begin(), edges_added[i].end(), e_cmp);
    //remove longest edges first (read from right to left)
    std::stable_sort(edges_removed[i].begin(), edges_removed[i].end(), e_cmp);
  }

  number_of_arrows += n; //count vertices additions
  edge_filtration.reserve(number_of_arrows);

//all vertices
  for(size_t i=0; i<n; ++i) {
    edge_filtration.emplace_back(i, i, std::numeric_limits<double>::infinity(), true);    
  }
  // edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add p_0
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    // edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(auto edg_it = edges_added[i].begin(); 
            edg_it != edges_added[i].end(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
    for(auto edg_it = edges_removed[i].rbegin(); 
             edg_it != edges_removed[i].rend(); ++edg_it) {
      edge_filtration.push_back(*edg_it);
    }
  }
}
#endif