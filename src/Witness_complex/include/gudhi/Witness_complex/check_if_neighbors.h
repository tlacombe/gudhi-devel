#ifndef CHECK_IF_NEIGHBORS_H_
#define CHECK_IF_NEIGHBORS_H_

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
    int stroke = 0; // absolute value = 0 if identical, 1 if one difference, 2 if more than one difference
    bool differ = false;
    auto svr1 = sc.simplex_vertex_range(simplex1).begin();
    auto svr2 = sc.simplex_vertex_range(simplex2).begin();
    auto svr1_end = sc.simplex_vertex_range(simplex1).end();
    auto svr2_end = sc.simplex_vertex_range(simplex2).end();    
    while (svr1 != svr1_end && svr2 != svr2_end && std::abs(stroke) != 2) {
      if (*svr1 == *svr2) {
        coface.push_back(*svr1);
        svr1++;
        svr2++;
      }
      else if (*svr1 > *svr2) {
        stroke++;
        coface.push_back(*(svr1++));
        differ = true;
      }
      else {
        stroke--;
        coface.push_back(*(svr2++));
        differ = true;
      }
    }
    if (std::abs(stroke) == 1 && svr1 != svr1_end) {
      coface.push_back(*svr1);
      return true;
    }
    if (std::abs(stroke) == 1 && svr2 != svr2_end) {
      coface.push_back(*svr2);
      return true;
    }
    return (std::abs(stroke) == 0 && differ);
  }

#endif
