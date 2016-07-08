#ifndef SAMPLING_RADIUS_H
#define SAMPLING_RADIUS_H

#include <vector>
#include <gudhi/distance_functions.h>
#include <CGAL/Delaunay_triangulation.h>


/* \brief Compute the sampling radius of a point set and a center realizing it.
 * \details The function takes as input a point range `points` and a predicate
 *  `pred` endowed with a boolean operator (), the purpose of which is to identify the balls that
 *  contribute to the sampling radius computation. 
 *  The output consists of a pair (sampling radius, center realizing it).
 */
template < typename K,
           typename Point_range,
           typename Predicate >
std::pair< typename K::FT, typename K::Point_d > sampling_radius(Point_range& points, Predicate pred)
{
  typedef typename K::FT       FT;
  typedef typename K::Point_d  Point_d;
  typedef typename K::Sphere_d Sphere_d;
  
  CGAL::Delaunay_triangulation<K> t(points[0].dimension());
  for (auto p: points)
    t.insert(p);
  FT epsilon2 = 0;
  Point_d final_center;
  Point_d control_point;
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it) {
    if (t.is_infinite(fc_it))
      continue;
    std::vector<Point_d> vertices;
    for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
      vertices.push_back((*fc_v_it)->point());
    Sphere_d cs(vertices.begin(), vertices.end());
    Point_d csc = cs.center();
    if (!pred(csc))
      continue;
    FT r2 = euclidean_distance(csc, *(vertices.begin()));
    if (epsilon2 < r2) {
      epsilon2 = r2;
      final_center = csc;         
      control_point = (*vertices.begin());
    }
  }
  return std::make_pair(sqrt(epsilon2), final_center);
}

/* \brief Compute the sampling radius of a point set and a center realizing it.
 * \details The function takes as input a point range `points`.
 *  The output consists of a pair (sampling radius, center realizing it).
 */
template < typename K >
std::pair< typename K::FT, typename K::Point_d > sampling_radius(std::vector<typename K::Point_d>& points)
{
  struct {
    bool operator ()(typename K::Point_d& csc)
    {
      bool in_ball = true; 
      for (auto xi = csc.cartesian_begin(); xi != csc.cartesian_end(); ++xi){
        	
		double dist=0;
		for(int r=0;r<csc.dimension();r++)
		{
			dist+=csc[r]*csc[r];		
		}
		if (dist>1) {
          	in_ball = false; break;
        	}}
      return in_ball;
    }
  } ball_pred;

  typedef typename K::Point_d Point_d;
  return sampling_radius<K>(points, ball_pred);
}

/* \brief Compute the sampling radius of a point set embedded on a flat torus and a center realizing it.
 * \details The function takes as input a point range `points` on a flat torus, representing the set $\mathbb{T}^d = (\mathbb{R}/[-1,1])^d$.
 *  The output consists of a pair (sampling radius, center realizing it).
 */
template < typename K >
std::pair< typename K::FT, typename K::Point_d > sampling_radius_torus(std::vector<typename K::Point_d>& points)
{
  struct {
    bool operator ()(typename K::Point_d& csc)
    {
      bool in_cube = true; 
      for (auto xi = csc.cartesian_begin(); xi != csc.cartesian_end(); ++xi)
        if (*xi > 1.0 || *xi < -1.0) {
          in_cube = false; break;
        }
      return in_cube;
    }
  } flat_torus_pred;

  return sampling_radius(points, flat_torus_pred);
}

#endif
