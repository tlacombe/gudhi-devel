#ifndef SAL_H
#define SAL_H

#include <list>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <limits>

namespace Gudhi {

typedef std::size_t Vertex;
typedef std::unordered_set<Vertex> Simplex;

std::list<Simplex> facets(const Simplex& sigma);
bool included(const Simplex& tau1, const Simplex& tau2);

class SAL {
    
public:
    void insert_max(const Simplex& sigma);

    void add(const Simplex& tau);
    void remove(const Simplex& tau);
    
    bool membership(const Simplex& tau) const;
    bool all_facets_inside(const Simplex& sigma);
    bool maximality(const Simplex& sigma) const;
    std::list<Simplex> max_cofaces(const Simplex& tau) const;

    //bool collapse(const Simplex& tau);
    Vertex contraction(const Vertex x, const Vertex y);

    std::size_t size() const;
    std::size_t num_vertices() const;
    
    typedef std::shared_ptr<Simplex> Simplex_ptr;
    struct Sptr_hash{ std::size_t operator()(const Simplex_ptr& s) const; };
    struct Sptr_equal{ std::size_t operator()(const Simplex_ptr& a, const Simplex_ptr& b) const; };
    typedef std::unordered_set<Simplex_ptr, Sptr_hash, Sptr_equal> Simplex_ptr_set;

private:
    void erase_max(const Simplex& sigma);
    Vertex best_index(const Simplex& tau) const;
    
    std::unordered_map<Vertex, Simplex_ptr_set> t0;
    bool max_empty_face; // Is the empty simplex a maximal face ?

};


/* sigma must not have neither face nor coface in the complex */
void SAL::insert_max(const Simplex& sigma){
    max_empty_face = (sigma.size()==0);
    Simplex_ptr sptr = std::make_shared<Simplex>(sigma);
    for(const Vertex& v : sigma){
        if(!t0.count(v)) t0.emplace(v, Simplex_ptr_set());
        t0.at(v).emplace(sptr);
    }
}

void SAL::add(const Simplex& tau){
    if(membership(tau)) return;
    bool facets_max = true;
    for(const Simplex& f : facets(tau))
        if(!maximality(f)) facets_max=false;
    if(facets_max)
        for(const Simplex& f : facets(tau))
            erase_max(f);
    else
        for(const Vertex& v : tau)
            //Copy constructor needed because the set is modified
            if(t0.count(v))  for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
                if(included(*sptr,tau)) erase_max(*sptr); // We erase all the maximal faces of tau
    insert_max(tau);
}

void SAL::remove(const Simplex& tau){
    if(tau.size()==0){
        t0.clear();
        max_empty_face = false;
    }
    else {
        const Vertex& v = best_index(tau);
        //Copy constructor needed because the set is modified
        if(t0.count(v)) for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
            if(included(tau, *sptr)){
                erase_max(*sptr);
                for(const Simplex& f : facets(tau))
                    if(!membership(f)) insert_max(f); // We add the facets which are new maximal simplices
            }
    }
}

bool SAL::membership(const Simplex& tau) const{
    if(t0.size()==0 && !max_empty_face) return false; //empty complex
    if(tau.size()==0) return true; //empty query simplex
    if(maximality(tau)) return true;
    const Vertex& v = best_index(tau);
    if(!t0.count(v))  return false;
    for(const Simplex_ptr& sptr : t0.at(v))
        if(included(tau, *sptr)) return true;
    return false;
}

bool SAL::all_facets_inside(const Simplex& sigma){
    Vertex v = best_index(sigma);
    if(!t0.count(v))  return false;
    Simplex f = sigma; f.erase(v);
    if(!membership(f)) return false;
    std::unordered_set<Vertex> facets_inside;
    for(const Simplex_ptr& sptr : t0.at(v))
        for(const Vertex& w : sigma){
            f = sigma; f.erase(w);
            if(included(f, *sptr)) facets_inside.insert(w);
        }
    return facets_inside.size() == sigma.size() - 1;
}


bool SAL::maximality(const Simplex& sigma) const{
    if(t0.size()==0 && !max_empty_face) return false; //empty complex
    if(sigma.size()==0) return max_empty_face;
    const Vertex& v =  best_index(sigma);
    if(!t0.count(v)) return false;
    return t0.at(v).count(std::make_shared<Simplex>(sigma));
}

std::list<Simplex> SAL::max_cofaces(const Simplex& tau) const{
    std::list<Simplex> cofaces;
    if(maximality(tau))
        cofaces.emplace_front(tau);
    else if(tau.size()==0){
        Simplex_ptr_set all_sptrs;
        for(const auto& kv : t0)
            for(const Simplex_ptr& sptr : kv.second) //kv.second is a Simplex_ptr_set
                all_sptrs.emplace(sptr);
        for(const Simplex_ptr& sptr : all_sptrs)
            cofaces.emplace_front(*sptr);
    }
    else {
        const Vertex& v = best_index(tau);
        if(t0.count(v)) for(const Simplex_ptr& sptr : t0.at(v))
            if(included(tau, *sptr)) cofaces.emplace_front(*sptr);
    }
    return cofaces;
}

/* Returns the remaining vertex */
Vertex SAL::contraction(const Vertex x, const Vertex y){
    if(!t0.count(x)) return y;
    if(!t0.count(y)) return x;
    int k, d;
    if(t0.at(x).size() > t0.at(y).size())
        k=x, d=y;
    else
        k=y, d=x;
    //Copy constructor needed because the set is modified
    for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(d))){
        Simplex sigma(*sptr);
        erase_max(sigma);
        sigma.erase(d);
        sigma.insert(k);
        add(sigma);
    }
    return k;
}

/* No facets added */
inline void SAL::erase_max(const Simplex& sigma){
    max_empty_face = false;
    Simplex_ptr sptr = std::make_shared<Simplex>(sigma);
    for(const Vertex& v : sigma){
        t0.at(v).erase(sptr);
        if(t0.at(v).size()==0) t0.erase(v);
    }
}


Vertex SAL::best_index(const Simplex& tau) const{
    std::size_t min = std::numeric_limits<size_t>::max();
    Vertex arg_min = -1;
    for(const Vertex& v : tau)
        if(!t0.count(v)) return v;
        else if(t0.at(v).size() < min)
            min = t0.at(v).size(), arg_min = v;
    return arg_min;
}

std::size_t SAL::size() const{
   Simplex s;
   return max_cofaces(s).size();
}

/* Return the number of vertices
 */
std::size_t SAL::num_vertices() const {
  return t0.size();
}

  
std::size_t SAL::Sptr_equal::operator()(const Simplex_ptr& s1, const Simplex_ptr& s2) const {
    if (s1->size() != s2->size()) return false;
    return included(*s1,*s2); //tests equality for same size simplices
}


std::size_t SAL::Sptr_hash::operator()(const Simplex_ptr& s) const {
    std::hash<double> h_f; //double hash is better than int hash
    size_t h = 0;
    for(const Vertex& v : *s)
        h += h_f(static_cast<double>(v));
    return h;
}

std::list<Simplex> facets(const Simplex& sigma){
    std::list<Simplex> facets;
    Simplex f(sigma);
    for(const Vertex& v : sigma){
        f.erase(v);
        facets.emplace_front(f);
        f.insert(v);
    }
    return facets;
}

bool included(const Simplex &tau1, const Simplex &tau2){
    for(const Vertex& v : tau1)
        if(!tau2.count(v)) return false;
    return true;
}
    
} //namespace Gudhi

#endif /* SAL_H */
