#pragma once
#include "Persistence_diagram.hpp"
#include <list>
#include <set>

// Neighbors_finder is an abstract data structure used to find if a query point has neighbors in the Neighbors_finder's persistence diagram.
// Persistence diagram's points have to be added manually using their index. A neighbor returned is automatically removed but we can also
// remove points manually using their index.
class Neighbors_finder{

protected:
    const Persistence_diagram& points;
    double r;

public:
    Neighbors_finder(const Persistence_diagram& points, double r);
    virtual ~Neighbors_finder() =0;
    virtual void add(int point_index) =0;
    virtual void remove(int point_index) =0;
    virtual bool contains(int point_index) const =0;
    virtual int pull_near(Diagram_point q) =0;
    virtual std::list<int>* pull_all_near(Diagram_point q);
};


// Naive_nf is a nave implementation of Neighbors_finder
class Naive_nf : public Neighbors_finder{

private:
    std::set<int> candidates;

public:
    Naive_nf(const Persistence_diagram& points, double r);
    void add(int point_index);
    void remove(int point_index);
    bool contains(int point_index) const;
    int pull_near(Diagram_point q);
};


// Simple_nf is the used Neighbors_finder's implementation
typedef Naive_nf Simple_nf;



Neighbors_finder::Neighbors_finder(const Persistence_diagram& points, double r) :
    points(points), r(r)
{}

inline Neighbors_finder::~Neighbors_finder()
{}

inline std::list<int>* Neighbors_finder::pull_all_near(Diagram_point q)
{
    std::list<int>* alln = new std::list<int>();
    int lastn = pull_near(q);
    while(lastn != null_point_index())
    {
        alln->emplace_back(lastn);
        lastn = pull_near(q);
    }
    return alln;
}

Naive_nf::Naive_nf(const Persistence_diagram & points, double r) :
    Neighbors_finder(points,r), candidates()
{}

void  Naive_nf::add(int point_index)
{
    candidates.emplace(point_index);
}

void Naive_nf::remove(int point_index)
{
    candidates.erase(point_index);
}

bool Naive_nf::contains(int point_index) const
{
    return (candidates.count(point_index) > 0);
}

int Naive_nf::pull_near(Diagram_point q)
{
    for(auto it = candidates.begin(); it != candidates.end(); ++it)
        if(distance(points.get_point(*it),q)<=r){
            int tmp = *it;
            candidates.erase(it);
            return tmp;
        }
    return null_point_index();
}
