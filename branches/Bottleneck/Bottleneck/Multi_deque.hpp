#pragma once
#include "Persistence_diagram.hpp"
#include <vector>
#include <list>

class Multi_deque{

private:
    int base_index;
    std::vector<int> prev;
    std::vector<int> next;
    std::vector<int> sentries_to_projections;
    std::vector<int> projections_to_sentries;

public:
    Multi_deque(int base_index, int projection_nb);
    void add_layer();
    void add(int point_index, int layer);
    void remove(int point_index);
    bool contains(int point_index, int layer) const;
    bool empty(int layer) const;
    int pull(int layer);
    std::list<int>* pull_all(int layer);
};



Multi_deque::Multi_deque(int base_index, int projections_nb) :
    base_index(base_index), prev(projections_nb, null_point_index()), next(projections_nb, null_point_index()),
    sentries_to_projections(1, null_point_index()), projections_to_sentries(projections_nb, -1)
{}

inline void Multi_deque::add_layer(){
    sentries_to_projections.emplace_back(null_point_index());
}

inline bool Multi_deque::contains(int point_index, int layer) const
{
    return projections_to_sentries.at(point_index - base_index)==layer;
}

inline bool Multi_deque::empty(int layer) const
{
    return sentries_to_projections.at(layer) == null_point_index();
}

inline void Multi_deque::remove(int point_index)
{
    point_index -= base_index;
    int layer = projections_to_sentries.at(point_index);
    if(point_index == sentries_to_projections.at(layer))
        sentries_to_projections[layer] = next.at(point_index);
    else{
        next[prev.at(point_index)] = next.at(point_index);
        if(next.at(point_index) != null_point_index())
            prev[next.at(point_index)] = prev.at(point_index);
    }
    projections_to_sentries.at(point_index) = -1;
}

inline void Multi_deque::add(int point_index, int layer)
{
    point_index -= base_index;
    projections_to_sentries[point_index] = layer;
    next[point_index] = sentries_to_projections.at(layer);
    if(sentries_to_projections.at(layer)!=null_point_index())
        prev[sentries_to_projections.at(layer)] = point_index;
    sentries_to_projections[layer] = point_index;
}

inline int Multi_deque::pull(int layer)
{
    int r = sentries_to_projections.at(layer);
    if(r!=null_point_index()){
        r+= base_index;
        remove(r);
    }
    return r;
}

inline std::list<int>* Multi_deque::pull_all(int layer)
{
    std::list<int>* all = new std::list<int>();
    for(int it = sentries_to_projections.at(layer); it!=null_point_index(); it = next.at(it))
        all->emplace_back(it);
    return all;
}
