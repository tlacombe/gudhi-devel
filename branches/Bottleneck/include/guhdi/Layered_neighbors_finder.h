#include "Neighbors_finder.h"

// Layered_neighbors_finder is a data structure used to find if a query point from U has neighbors in V in a given layer of the
// layered persistence diagrams graph. V's points have to be added manually using their index. A neighbor returned is automatically
// removed.
class Layered_neighbors_finder{

public:
    Layered_neighbors_finder(const Persistence_diagrams_graph& g, double r);
    void add(int v_point_index, int layer);
    int pull_near(int u_point_index, int layer);
    int layers_number() const;

private:
    const Persistence_diagrams_graph& g;
    const double r;
    std::vector<Neighbors_finder> neighbors_finder;
};



Layered_neighbors_finder::Layered_neighbors_finder(const Persistence_diagrams_graph& g, double r) :
    g(g), r(r), neighbors_finder()
{}

inline void Layered_neighbors_finder::add(int v_point_index, int layer)
{
    for(int l = neighbors_finder.size(); l<layer; l++)
        neighbors_finder.emplace_back(Neighbors_finder(g,r));
    neighbors_finder.at(layer).add(v_point_index);
}

inline int Layered_neighbors_finder::pull_near(int u_point_index, int layer)
{
    return neighbors_finder.at(layer).pull_near(u_point_index);
}

inline int Layered_neighbors_finder::layers_number() const
{
    return neighbors_finder.size();
}
