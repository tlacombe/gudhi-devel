#pragma once
#include "Neighbors_finder.hpp"
#include <vector>

// Multi_nf is a data structure used to find if a query point has neighbors in a layered persistence diagram.
// Layer's points have to be added manually using their index. A neighbor returned is automatically removed but we can also
// remove points manually using their index.
class Multi_nf{

private:
    const Persistence_diagram& points;
    // Multi_nf has a Neighbors_finder for each layer.
    std::vector<Simple_nf> layers_to_nf;
    std::vector<int> points_to_layers;
    double r;

public:
    Multi_nf(const Persistence_diagram & p, double r);
    void add_layer();
    void add(int point_index, int layer);
    void remove(int point_index);
    bool contains(int point_index, int layer) const;
    int pull_near(Diagram_point q, int layer);
    std::list<int>* pull_all_near(Diagram_point q, int layer);
};



Multi_nf::Multi_nf(const Persistence_diagram &p, double r) :
    points(p), layers_to_nf(1, Simple_nf(p, r)), points_to_layers(p.size(),-1), r(r){}

inline void Multi_nf::add_layer()
{
    layers_to_nf.emplace_back(Simple_nf(points, r));
}

inline void Multi_nf::add(int point_index, int layer)
{
    layers_to_nf.at(layer).add(point_index);
    points_to_layers[point_index] = layer;
}

inline void Multi_nf::remove(int point_index)
{
    int l = points_to_layers.at(point_index);
    layers_to_nf.at(l).remove(point_index);
    points_to_layers[point_index] = -1;
}

inline bool Multi_nf::contains(int point_index, int layer) const
{
    return points_to_layers.at(point_index) == layer;
}

inline int Multi_nf::pull_near(Diagram_point q, int layer)
{
    int r = layers_to_nf.at(layer).pull_near(q);
    if(r != null_point_index())
        points_to_layers[r] = -1;
    return r;
}

inline std::list<int>* Multi_nf::pull_all_near(Diagram_point q, int layer)
{
    std::list<int>* alln = layers_to_nf.at(layer).pull_all_near(q);
    for(auto it = alln->cbegin(); it != alln->cend(); ++it)
        points_to_layers[*it] = -1;
    return alln;
}

