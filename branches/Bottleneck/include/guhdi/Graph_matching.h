#include "Layered_neighbors_finder.h"
#include <deque>

//template<typename Persistence_diagram1, typename Persistence_diagram2>
//double bottleneck_distance(Persistence_diagram1& diag1, Persistence_diagram2& diag2, double e = 0.);


class Graph_matching{

public:
    Graph_matching(const Persistence_diagrams_graph& g);
    Graph_matching& operator=(const Graph_matching& m);
    bool perfect() const;
    bool multi_augment();
    void set_r(double r);

private:
    const Persistence_diagrams_graph& g;
    double r;
    std::vector<int> v_to_u;
    std::list<int> unmatched_in_u;


    Layered_neighbors_finder* layering() const;
    bool augment(Layered_neighbors_finder* layered_nf, int start, int max_depth);
    void update(std::deque<int>& path);
};

Graph_matching::Graph_matching(const Persistence_diagrams_graph& g)
    : g(g), r(0), v_to_u(g.size()), unmatched_in_u()
{
    for(int u_point_index=0; u_point_index<g.size(); ++u_point_index)
        unmatched_in_u.emplace_back(u_point_index);
}

Graph_matching& Graph_matching::operator=(const Graph_matching& m)
{
    r = m.r;
    v_to_u = m.v_to_u;
    unmatched_in_u = m.unmatched_in_u;
    return *this;
}

inline bool Graph_matching::perfect() const
{
    return unmatched_in_u.empty();
}

inline bool Graph_matching::multi_augment()
{
    if(perfect())
        return false;
    Layered_neighbors_finder* layered_nf = layering();
    double rn = sqrt(g.size());
    int ln = layered_nf->layers_number();
    // verification of a necessary criterion
    if((unmatched_in_u.size() > rn && ln*2. >= rn) || ln==0)
        return false;
    std::list<int>* tries = new std::list<int>(unmatched_in_u);
    for(auto it = tries->cbegin(); it != tries->cend(); it++)
        augment(layered_nf, *it, ln);
    delete tries;
    delete layered_nf;
    // Since layers_number > 0, we know that at least 1 augment has returned true
    return true;
}

inline void Graph_matching::set_r(double r){
    this->r = r;
}

Layered_neighbors_finder* Graph_matching::layering() const
{
    bool end = false;
    int layer = 0;
    std::list<int> u_vertices(unmatched_in_u);
    std::list<int> v_vertices;
    Neighbors_finder nf(g,r);
    Layered_neighbors_finder* layered_nf = new Layered_neighbors_finder(g,r);
    for(int v_point_index=0; v_point_index<g.size(); ++v_point_index)
        nf.add(v_point_index);
    while(!u_vertices.empty()){
        for(auto it = u_vertices.cbegin(); it != u_vertices.cend(); ++it){
            std::list<int>* u_succ = nf.pull_all_near(*it);
            for(auto it = u_succ->cbegin(); it != u_succ->cend(); ++it){
                layered_nf->add(*it, layer);
                v_vertices.emplace_back(*it);
            }
            delete u_succ;
        }
        u_vertices.clear();
        for(auto it = v_vertices.cbegin(); it != v_vertices.cend(); it++){
            if(v_to_u.at(*it)==null_point_index())
                end = true;
            else
                u_vertices.emplace_back(v_to_u.at(*it));
        }
        if(end)
            return layered_nf;
        v_vertices.clear();
        layer++;
    }
    return layered_nf;
}

bool Graph_matching::augment(Layered_neighbors_finder *layered_nf, int start, int max_depth)
{
    std::deque<int> path;
    path.emplace_back(start);
    do{
        if((int) path.size() > max_depth*2 +1){
            path.pop_back();
            path.pop_back();
        }
        if(path.empty())
            return false;
        path.emplace_back(layered_nf->pull_near(path.back(),path.size()/2));
        while(path.back()==null_point_index()){
            path.pop_back();
            path.pop_back();
            if(path.empty())
                return false;
            path.pop_back();
            path.emplace_back(layered_nf->pull_near(path.back(),path.size()/2));
        }
        path.emplace_back(v_to_u.at(path.back()));
    }
    while(path.back()!=null_point_index());
    path.pop_back();
    update(path);
    return true;
}

void Graph_matching::update(std::deque<int>& path)
{
    unmatched_in_u.remove(path.front());
    for(auto it = path.cbegin(); it != path.cend(); ++it){
        int tmp = *it;
        ++it;
        v_to_u[*it] = tmp;
    }
}

template<typename Persistence_diagram1, typename Persistence_diagram2>
double bottleneck_distance(Persistence_diagram1& diag1, Persistence_diagram2& diag2, double e)
{
    Persistence_diagrams_graph g(diag1, diag2, e);
    std::vector<double>* sd = g.sorted_distances();
    int idmin = 0;
    int idmax = sd->size()-1;
    double alpha = pow(sd->size(), 0.25);
    Graph_matching m(g);
    Graph_matching biggest_unperfect =m;
    while(idmin != idmax){
        int pas = (int) ((idmax-idmin)/alpha);
        m.set_r(sd->at(idmin + pas));
        while(m.multi_augment());
        if(m.perfect()){
            idmax = idmin + pas;
            m = biggest_unperfect;
        }
        else{
            biggest_unperfect = m;
            idmin = idmin + pas + 1;
        }
    }
    double b = sd->at(idmin);
    delete sd;
    return b;
}
