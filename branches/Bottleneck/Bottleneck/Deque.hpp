#pragma once
#include "Multi_deque.hpp"

class Deque{

private:
    Multi_deque md;

public:
    Deque(int base_index, int projection_nb);
    void add(int point_index);
    void remove(int point_index);
    bool contains(int point_index);
    bool empty();
    int pull();
    std::list<int>* pull_all();
};



Deque::Deque(int base_index, int projection_nb):
    md(base_index, projection_nb)
{
md.add_layer();
}

inline void Deque::add(int point_index)
{
    md.add(point_index,0);
}

inline void Deque::remove(int point_index)
{
    md.remove(point_index);
}

inline bool Deque::contains(int point_index)
{
    return md.contains(point_index,0);
}

inline bool Deque::empty()
{
    return md.empty(0);
}

inline int Deque::pull()
{
    return md.pull(0);
}

inline std::list<int>* Deque::pull_all()
{
    return md.pull_all(0);
}
