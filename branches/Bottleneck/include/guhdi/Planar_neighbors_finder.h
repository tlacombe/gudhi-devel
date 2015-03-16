#include "Persistence_diagrams_graph.h"
#include <list>

#include <iostream>

// Planar_neighbors_finder is a data structure used to find if a query point from U has planar neighbors in V with the planar distance.
// V's points have to be added manually using their index. A neighbor returned is automatically removed but we can also
// remove points manually using their index.
class Abstract_planar_neighbors_finder{

public:
    Abstract_planar_neighbors_finder(const Persistence_diagrams_graph& g, double r);
    virtual ~Abstract_planar_neighbors_finder() =0;
    virtual void add(int v_point_index) =0;
    virtual void remove(int v_point_index) =0;
    virtual bool contains(int v_point_index) const =0;
    virtual int pull_near(int u_point_index) =0;
    virtual std::list<int>* pull_all_near(int u_point_index);

protected:
    const Persistence_diagrams_graph& g;
    const double r;
};


// Naive_pnf is a nave implementation of Abstract_planar_neighbors_finder
class Naive_pnf : public Abstract_planar_neighbors_finder{

public:
    Naive_pnf(const Persistence_diagrams_graph& g, double r);
    void add(int v_point_index);
    void remove(int v_point_index);
    bool contains(int v_point_index) const;
    int pull_near(int u_point_index);

private:
    std::set<int> candidates;
};


// Planar_neighbors_finder is the used Abstract_planar_neighbors_finder's implementation
typedef Naive_pnf Planar_neighbors_finder;



Abstract_planar_neighbors_finder::Abstract_planar_neighbors_finder(const Persistence_diagrams_graph& g, double r) :
    g(g), r(r)
{}

inline Abstract_planar_neighbors_finder::~Abstract_planar_neighbors_finder()
{}

inline std::list<int>* Abstract_planar_neighbors_finder::pull_all_near(int u_point_index)
{
    std::list<int>* all_pull = new std::list<int>();
    int last_pull = pull_near(u_point_index);
    while(last_pull != null_point_index())
    {
        all_pull->emplace_back(last_pull);
        last_pull = pull_near(u_point_index);
    }
    return all_pull;
}

Naive_pnf::Naive_pnf(const Persistence_diagrams_graph& g, double r)  :
    Abstract_planar_neighbors_finder(g,r), candidates()
{}

inline void  Naive_pnf::add(int v_point_index)
{
    candidates.emplace(v_point_index);
}

inline void Naive_pnf::remove(int v_point_index)
{
    candidates.erase(v_point_index);
}

inline bool Naive_pnf::contains(int v_point_index) const
{
    return (candidates.count(v_point_index) > 0);
}

inline int Naive_pnf::pull_near(int u_point_index)
{
    for(auto it = candidates.begin(); it != candidates.end(); ++it)
        if(g.distance(u_point_index,*it)<=r){
            int tmp = *it;
            candidates.erase(it);
            return tmp;
        }
    return null_point_index();
}
