#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

namespace subsampl = Gudhi::subsampling;
namespace tc = Gudhi::tangential_complex;

typedef CGAL::Epick_d<CGAL::Dimension_tag<3>>                   Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_d                                         Point;
typedef Kernel::Vector_d                                        Vector;
typedef tc::Tangential_complex<
  Kernel, CGAL::Dimension_tag<2>,
  CGAL::Parallel_tag>                                           TC;

int main (void)
{
  const int INTRINSIC_DIM = 2;
  const int AMBIENT_DIM = 3;
  const int NUM_POINTS = 50;

  Kernel k;

  // Generate points on a 2-sphere
  CGAL::Random_points_on_sphere_d<Point> generator(AMBIENT_DIM, 3.);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i < NUM_POINTS; ++i)
    points.push_back(*generator++);

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, 0.05, k);
  tc.compute_tangential_complex();

  // Try to fix inconsistencies
  unsigned int num_perturb_steps;
  std::size_t initial_num_inconsistent_local_tr;
  std::size_t best_num_inconsistent_local_tr;
  std::size_t final_num_inconsistent_local_tr;
  auto perturb_ret = tc.fix_inconsistencies_using_perturbation(
    num_perturb_steps, initial_num_inconsistent_local_tr,
    best_num_inconsistent_local_tr, final_num_inconsistent_local_tr,
    10); // give it 10 seconds to succeed

  // Export the TC into a Simplex_tree
  Gudhi::Simplex_tree<> stree;
  tc.export_complex(stree);

  return 0;
}
