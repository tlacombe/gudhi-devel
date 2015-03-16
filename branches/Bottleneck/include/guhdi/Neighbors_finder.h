#include "Planar_neighbors_finder.h"
#include <unordered_set>

// Neighbors_finder is a data structure used to find if a query point from U has neighbors in V in the persistence diagrams graph.
// V's points have to be added manually using their index. A neighbor returned is automatically removed.
class Neighbors_finder{

public:
    Neighbors_finder(const Persistence_diagrams_graph& g, double r);
    void add(int v_point_index);
    int pull_near(int u_point_index);
    std::list<int>* pull_all_near(int u_point_index);

private:
    const Persistence_diagrams_graph& g;
    const double r;
    Planar_neighbors_finder planar_neighbors_f;
    std::unordered_set<int> projections_f;
};



Neighbors_finder::Neighbors_finder(const Persistence_diagrams_graph& g, double r) :
    g(g), r(r), planar_neighbors_f(g,r), projections_f()
{}

inline void Neighbors_finder::add(int v_point_index)
{
    if(g.on_the_v_diagonal(v_point_index))
        projections_f.emplace(v_point_index);
    else
        planar_neighbors_f.add(v_point_index);
}

inline int Neighbors_finder::pull_near(int u_point_index)
{
    int v_challenger = g.corresponding_point_in_v(u_point_index);
    if(planar_neighbors_f.contains(v_challenger) && g.distance(u_point_index, v_challenger) < r){
         planar_neighbors_f.remove(v_challenger);
         return v_challenger;
    }
    if(g.on_the_u_diagonal(u_point_index)){
        auto it = projections_f.cbegin();
        if(it != projections_f.cend()){
            int tmp = *it;
            projections_f.erase(it);
            return tmp;
        }
    }
    else
        return planar_neighbors_f.pull_near(u_point_index);
    return null_point_index();
}

inline std::list<int>* Neighbors_finder::pull_all_near(int u_point_index)
{
    std::list<int>* all_pull = planar_neighbors_f.pull_all_near(u_point_index);
    int last_pull = pull_near(u_point_index);
    while(last_pull != null_point_index())
    {
        all_pull->emplace_back(last_pull);
        last_pull = pull_near(u_point_index);
    }
    return all_pull;
}
