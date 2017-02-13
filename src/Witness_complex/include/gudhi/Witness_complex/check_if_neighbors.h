#ifndef CHECK_IF_NEIGHBORS_H_
#define CHECK_IF_NEIGHBORS_H_

/* \brief Check if two k-dimensional simplices have
   a common facet in a simplicial complex.
 */
template < typename SimplicialComplexForWitness,
           typename Simplex,
           typename VertexVector >
bool check_if_neighbors(SimplicialComplexForWitness& sc,
                        Simplex& simplex1,
                        Simplex& simplex2,
                        VertexVector& coface
                        )
  {
    int stroke = 0; // 0 is identical, 1 is one difference, 2 is more than one difference
    auto svr1 = sc.simplex_vertex_range(simplex1).begin();
    auto svr2 = sc.simplex_vertex_range(simplex2).begin();
    auto svr1_end = sc.simplex_vertex_range(simplex1).end();
    auto svr2_end = sc.simplex_vertex_range(simplex2).end();    
    while (svr1 != svr1_end && svr2 != svr2_end && stroke != 2) {
      if (*svr1 == *svr2) {
        coface.push_back(*svr1);
        svr1++;
        svr2++;
      }
      else {
        stroke++;
        if (*svr1 > *svr2)
          coface.push_back(*(svr1++));
        else
          coface.push_back(*(svr2++));
      }
    }
    if (svr1 != svr1_end)
      coface.push_back(*svr1);
    if (svr2 != svr2_end)
      coface.push_back(*svr2);
    return (stroke == 1);
  }

#endif
