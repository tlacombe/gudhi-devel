#pragma once
#include "Multi_nf.hpp"
#include "Multi_deque.hpp"
#include "Deque.hpp"
#include <vector>
#include <cmath>
#include <deque>

class Graph_matching{

private:
    const Persistence_diagram & x0;
    const Persistence_diagram & y0;
    int n;
    double r;
    std::vector< int > v_to_u;
    Deque unmatched_in_u;

    int size();
    int layering(Multi_nf * mpnf, Multi_deque * md);
    bool augment(Multi_nf * mpnf, Multi_deque * md, int start_i, int d);
    int get_succ_in_v(int u_i, Simple_nf* mpnf, Deque* md);
    int get_succ_in_v(int u_i, Multi_nf * mpnf, Multi_deque * md, int d);
    void update(std::deque< int > & path);
    static std::vector<double>* sorted_distances(const Persistence_diagram & p1, const Persistence_diagram & p2);
    static double bottleneck_distance(const Persistence_diagram & p1, const Persistence_diagram & p2);

public:
    Graph_matching(const Persistence_diagram & x0, const Persistence_diagram & y0);
    Graph_matching(const Graph_matching & m);
    Graph_matching& operator=(const Graph_matching& m);
    void set_r(double r);
    bool multi_augment();
    bool perfect();

    template< typename Container_type1, typename Container_type2 >
    static double bottleneck_distance(Container_type1 & c1, Container_type2 & c2, double e = 0.);
};



    int Graph_matching::size()
    {
        int r = 0;
        for(int i=0; i<n; i++)
            if(v_to_u.at(i)!=null_point_index())
                r++;
        return r;
    }

    int Graph_matching::layering(Multi_nf* mpnf, Multi_deque* md)
    {
        bool end = false;
        int d = -1;
        std::list<int>* u_vertices(unmatched_in_u.pull_all());
        std::list<int> v_vertices;
        Simple_nf* mpnf0 = new Simple_nf(y0,r);
        Deque* md0 = new Deque(y0.size(),x0.size());
        for(int p_i = 0; p_i < y0.size(); p_i++)
            mpnf0->add(p_i);
        for(int p_i = y0.size(); p_i < x0.size() + y0.size(); p_i++)
            md0->add(p_i);
        while(!u_vertices->empty()){
            d++;
            for(auto it = u_vertices->cbegin(); it != u_vertices->cend(); it++)
                if(*it < x0.size()){
                    Diagram_point u_point = x0.get_point(*it);
                    int projection_i = *it + y0.size();
                    if(u_point.second - u_point.first <= 2*r && md0->contains(projection_i)){
                        md0->remove(projection_i);
                        md->add(projection_i, d);
                        v_vertices.emplace_back(projection_i);
                    }
                          std::list<int>* alln = mpnf0->pull_all_near(u_point);
                    for(auto it = alln->cbegin(); it != alln->cend(); ++it){
                        mpnf->add(*it, d);
                        v_vertices.emplace_back(*it);
                    }
                    delete alln;
                }
                else{
                    int next = get_succ_in_v(*it, mpnf0, md0);
                    while(next!=null_point_index()){
                        v_vertices.emplace_back(next);
                        if(next < y0.size())
                            mpnf->add(next, d);
                        else
                            md->add(next, d);
                        next = get_succ_in_v(*it, mpnf0, md0);
                    }
                }
            u_vertices->clear();
            for(auto it = v_vertices.cbegin(); it != v_vertices.cend(); it++){
                if(v_to_u.at(*it)==null_point_index())
                    end = true;
                else
                    u_vertices->emplace_back(v_to_u.at(*it));
            }
            if(end){
                delete mpnf0;
                delete md0;
                return d;
            }
            v_vertices.clear();
            mpnf->add_layer();
            md->add_layer();
        }
        delete mpnf0;
        delete md0;
        return -1;
    }

    bool Graph_matching::augment(Multi_nf* mpnf, Multi_deque* md, int start_i, int d)
    {
        std::deque< int > path;
        path.emplace_back(start_i);
        do{
            if((int) path.size() > d*2 +1){
                path.pop_back();
                path.pop_back();
            }
            if(path.empty())
                return false;
            path.emplace_back(get_succ_in_v(path.back(), mpnf, md, path.size()/2));
            while(path.back()==null_point_index()){
                path.pop_back();
                path.pop_back();
                if(path.empty())
                    return false;
                path.pop_back();
                path.emplace_back(get_succ_in_v(path.back(), mpnf, md, path.size()/2));
            }
            path.emplace_back(v_to_u.at(path.back()));
        }
        while(path.back()!=null_point_index());
        path.pop_back();
        update(path);
        return true;
    }

    int Graph_matching::get_succ_in_v(int u_i, Simple_nf* mpnf, Deque* md)
    {
        if(u_i < x0.size()){
            Diagram_point u_point = x0.get_point(u_i);
            int projection_i = u_i + y0.size();
            if(u_point.second - u_point.first <= 2*r && md->contains(projection_i)){
                md->remove(projection_i);
                return projection_i;
            }
            return mpnf->pull_near(u_point);
        }
        int op = md->pull();
        if(op != null_point_index())
            return op;
        int projector_i = u_i - x0.size();
        Diagram_point projector = y0.get_point(projector_i);
        if(projector.second - projector.first <= 2*r && mpnf->contains(projector_i)){
            mpnf->remove(projector_i);
            return projector_i;
        }
        return null_point_index();
    }


    int Graph_matching::get_succ_in_v(int u_i, Multi_nf * mpnf, Multi_deque * md, int d)
    {
        if(u_i < x0.size()){
            Diagram_point u_point = x0.get_point(u_i);
            int projection_i = u_i + y0.size();
            if(u_point.second - u_point.first <= 2*r && md->contains(projection_i, d)){
                md->remove(projection_i);
                return projection_i;
            }
            return mpnf->pull_near(u_point, d);
        }
        int op = md->pull(d);
        if(op != null_point_index())
            return op;
        int projector_i = u_i - x0.size();
        Diagram_point projector = y0.get_point(projector_i);
        if(projector.second - projector.first <= 2*r && mpnf->contains(projector_i, d)){
            mpnf->remove(projector_i);
            return projector_i;
        }
        return null_point_index();
    }

    void Graph_matching::update(std::deque< int > & path)
    {
        unmatched_in_u.remove(path.front());
        for(auto it = path.cbegin(); it != path.cend(); ++it){
            int tmp = *it;
            ++it;
            v_to_u[*it] = tmp;
        }
    }

    std::vector<double>* Graph_matching::sorted_distances(const Persistence_diagram & p1, const Persistence_diagram & p2)
    {
        double max = 0;
        Diagram_point p;
        for(int i = 0; i< p1.size(); i++){
            p = p1.get_point(i);
            if(p.second - p.first>max)
                max = p.second - p.first;
        }
        for(int i = 0; i< p2.size(); i++){
            p = p2.get_point(i);
            if(p.second - p.first>max)
                max = p.second - p.first;
        }
        std::vector<double>* sorted_distances = new std::vector<double>();
        for(double d = 0; d<= max/2; d += 1./2.)
            sorted_distances->emplace_back(d);
        return sorted_distances;
    }

    double Graph_matching::bottleneck_distance(const Persistence_diagram & p1, const Persistence_diagram & p2)
    {
        std::vector<double>* sd = sorted_distances(p1,p2);
        int idmin = 0;
        int idmax = sd->size()-1;
        double alpha = pow(sd->size(), 0.25);
        Graph_matching biggest_unperfect(p1,p2);
        Graph_matching m(biggest_unperfect);
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

    Graph_matching::Graph_matching(const Persistence_diagram & x0, const Persistence_diagram & y0) :
        x0(x0), y0(y0), n(x0.size() + y0.size()), r(), v_to_u(n, null_point_index()), unmatched_in_u(0,n)
    {
        for(int i=0; i<n; i++)
            unmatched_in_u.add(i);
    }

    Graph_matching::Graph_matching(const Graph_matching & m):
        x0(m.x0), y0(m.y0), n(m.n), r(m.r), v_to_u(m.v_to_u), unmatched_in_u(m.unmatched_in_u)
    {}

    Graph_matching& Graph_matching::operator=(const Graph_matching& m)
    {
        n = m.n;
        r = m.r;
        v_to_u = m.v_to_u;
        unmatched_in_u = m.unmatched_in_u;
        return *this;
    }

    void Graph_matching::set_r(double r)
    {
        this->r=r;
    }

    bool Graph_matching::multi_augment()
    {
        if(perfect())
            return false;
        Multi_nf * mpnf = new Multi_nf(y0,r);
        Multi_deque * md = new Multi_deque(y0.size(),x0.size());
        int d = layering(mpnf, md);
        double rn = sqrt(n);
        if((size()< n - rn && d*2 >= rn) || d<0)
            return false;
        std::list<int>* tries = unmatched_in_u.pull_all();
        for(auto it = tries->cbegin(); it != tries->cend(); it++)
            augment(mpnf, md, *it, d);
        delete tries;
        delete mpnf;
        delete md;
        return true;
    }

    bool Graph_matching::perfect()
    {
        return unmatched_in_u.empty();
    }


    template< typename Container_type1, typename Container_type2 >
    double Graph_matching::bottleneck_distance(Container_type1 & c1, Container_type2 & c2, double e){
        const Persistence_diagram & p1 = Persistence_diagram(c1, e);
        const Persistence_diagram & p2 = Persistence_diagram(c2, e);
        return p1.size() > p2.size() ? bottleneck_distance(p1, p2) : bottleneck_distance(p2, p1) ;
    }
