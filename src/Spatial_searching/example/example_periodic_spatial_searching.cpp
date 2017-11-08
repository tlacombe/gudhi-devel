#include <gudhi/Periodic_kd_tree_search.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <vector>
#include <string>
#include <sstream>

namespace gss = Gudhi::spatial_searching;

const int DIM = 2;

template <typename Point>
std::string print_point(Point const& p)
{
  std::stringstream sstr;
  sstr << "(";
  for (int i = 0; i < DIM - 1; ++i)
    sstr << p[i] << ", ";
  sstr << p[DIM - 1] << ")";
  return sstr.str();
}

int main(void) {
  typedef CGAL::Epick_d<CGAL::Dimension_tag<2> > K;
  typedef typename K::Point_d Point;
  typedef std::vector<Point> Points;

  typedef gss::Periodic_kd_tree_search<K, Points> Points_ds;

  Points points;
  points.push_back(Point(0., 0.));
  points.push_back(Point(0., 0.1));
  points.push_back(Point(0.9, 0.9));

  Points_ds points_ds(points, K::Iso_box_d(Point(0., 0.), Point(1., 1.) ));

  Points queries;
  queries.push_back(Point(0.99, 0.99));

  for (auto const& q : queries)
  {
    std::cout << "Query: " << print_point(q) << "\n";

    // 10-nearest neighbors
    std::cout << "10 nearest neighbors:\n";
    auto knn_range = points_ds.query_k_nearest_neighbors(q, 10, true);
    for (auto const& nghb : knn_range)
    {
      std::cout << "  " << nghb.first
        << " - " << print_point(points[nghb.first])
        << " - sq. dist. = " << nghb.second << "\n";
    }

    // 10-farthest neighbor query
    std::cout << "10 farthest neighbors:\n";
    auto kfn_range = points_ds.query_k_farthest_neighbors(q, 10, true);
    for (auto const& nghb : kfn_range)
    {
      std::cout << "  " << nghb.first
        << " - " << print_point(points[nghb.first])
        << " - sq. dist. = " << nghb.second << "\n";
    }

    // All-near-neighbors
    // TODO: see Periodic_kd_tree_search.h
    /*std::cout << "All-near-neighbors (radius = 0.2):\n";
    std::vector<std::size_t> rs_result;
    points_ds.near_search(q, 0.2, std::back_inserter(rs_result));
    K k;
    for (auto const& p_idx : rs_result)
      std::cout << "  " << p_idx << " (sq. dist. = " << k.squared_distance_d_object()(points[p_idx], q) << ")\n";*/
  }

  return 0;
}
