#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <gudhi/Landmark_choice_by_random_point.h>
#include <vector>
#include <iterator>

int main() {
  typedef int Point_d;
  std::vector<Point_d> vect = {0,1,2,3,4,5,6,7,8,9,10};
  std::vector<Point_d> landmarks;
  Gudhi::landmark_choice_by_random_point(vect, 5, std::back_inserter(landmarks));
  std::cout << "landmark vector contains: ";
  for (auto l: landmarks)
    std::cout << l << " ";
  std::cout << "\n";
}
