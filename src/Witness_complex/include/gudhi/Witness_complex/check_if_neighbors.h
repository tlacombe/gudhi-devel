#ifndef CHECK_IF_NEIGHBORS_H_
#define CHECK_IF_NEIGHBORS_H_

#include <set>

/* \brief Check if two k-dimensional simplices have
   a common facet in a simplicial complex.
 */
template < typename SimplicialComplexForWitness,
           typename Simplex,
           typename VertexVector >
bool check_if_neighbors(SimplicialComplexForWitness& sc,
                        const Simplex& simplex1,
                        const Simplex& simplex2,
                        VertexVector& coface
                        )
  {
    int stroke = 0;
    int diff_count = 0;
    auto svr1 = sc.simplex_vertex_range(simplex1).begin();
    auto svr2 = sc.simplex_vertex_range(simplex2).begin();
    auto svr1_end = sc.simplex_vertex_range(simplex1).end();
    auto svr2_end = sc.simplex_vertex_range(simplex2).end();    
    while (svr1 != svr1_end && svr2 != svr2_end && std::abs(stroke) != 2 && diff_count != 2) {
      if (*svr1 == *svr2) {
        coface.push_back(*svr1);
        svr1++;
        svr2++;
      }
      else if (*svr1 > *svr2) {
        coface.push_back(*(svr1++));
        if (++stroke != 0)
          diff_count++;
      }
      else {
        coface.push_back(*(svr2++));
        if (--stroke != 0)
          diff_count++;
      }
    }
    if (std::abs(stroke) == 1  && diff_count != 2 && svr1 != svr1_end) {
      coface.push_back(*svr1);
      return true;
    }
    if (std::abs(stroke) == 1 && diff_count != 2 && svr2 != svr2_end) {
      coface.push_back(*svr2);
      return true;
    }
    return (std::abs(stroke) == 0 && diff_count == 1);
  }

// template < typename SimplicialComplexForWitness,
//            typename Simplex,
//            typename VertexVector >
// bool check_if_neighbors(SimplicialComplexForWitness& sc,
//                         const Simplex& simplex1,
//                         const Simplex& simplex2,
//                         VertexVector& coface
//                         )
//   {
//     std::size_t k = sc.dimension(simplex1);
//     std::set<typename SimplicialComplexForWitness::Vertex_handle> vertex_set;
//     for (auto v: sc.simplex_vertex_range(simplex1))
//       vertex_set.insert(v);
//     for (auto v: sc.simplex_vertex_range(simplex2))
//       vertex_set.insert(v);
//     if (vertex_set.size() == k + 2) {
//       coface = VertexVector(vertex_set.begin(), vertex_set.end());
//       return true;
//     }
//     return false;
//   }


#endif
