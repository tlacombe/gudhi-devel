#ifndef SAL_H
#define SAL_H

#include <gudhi/Simplex.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <limits>

#define default_vertex 0
#define default_filtration 0

namespace Gudhi {

typedef Simplex::Vertex Vertex;

class SAL {

public:
    typedef std::shared_ptr<Simplex> Simplex_ptr;
    typedef Simplex::Filtration_t Filtration_value;

    struct Sptr_hash{ std::size_t operator()(const Simplex_ptr& s) const; };
    struct Sptr_equal{ std::size_t operator()(const Simplex_ptr& a, const Simplex_ptr& b) const; };
    typedef std::unordered_set<Simplex_ptr, Sptr_hash, Sptr_equal> Simplex_ptr_set;

    template <typename Input_vertex_range>
    void insert_critical_simplex(const Input_vertex_range &vertex_range, Filtration_value f = default_filtration);

    template <typename Input_vertex_range>
    void insert_simplex(const Input_vertex_range &vertex_range, Filtration_value f = default_filtration);

    template <typename Input_vertex_range>
    void remove_simplex(const Input_vertex_range &vertex_range);
    
    template <typename Input_vertex_range>
    Simplex_ptr find(const Input_vertex_range &vertex_range) const;

    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range) const;

    template <typename Input_vertex_range>
    bool criticality(const Input_vertex_range &vertex_range) const;

    template <typename Input_vertex_range>
    Simplex_ptr_set critical_cofaces(const Input_vertex_range &vertex_range) const;

    template <typename Input_vertex_range>
    Filtration_value get_filtration(const Input_vertex_range &vertex_range) const;

    Vertex contraction(const Vertex x, const Vertex y);

    std::size_t num_simplices() const;
    std::size_t num_vertices() const;
    

    // Simplex_tree interface compatibility
    typedef Simplex_ptr Simplex_handle;
    typedef Vertex Vertex_handle;
    typedef Simplex_ptr_set Complex_simplex_range;
    typedef Simplex Simplex_vertex_range; //Iterator over vertices of a simplex. More...
    Simplex_handle null_simplex() const;
    Simplex_vertex_range simplex_vertex_range (Simplex_handle const &simplex) const;

    // For Siargey
    template <typename Input_vertex_range>
    bool all_facets_inside(const Input_vertex_range &vertex_range) const;
    Simplex_ptr_set candidates() const;

protected:
    void erase_critical(const Simplex_ptr& sptr);

    template <typename Input_vertex_range>
    Vertex best_index(const Input_vertex_range &vertex_range) const;
    
    std::unordered_map<Vertex, Simplex_ptr_set> t0;
};

typedef SAL::Simplex_ptr Simplex_ptr;
typedef SAL::Simplex_ptr_set Simplex_ptr_set;

template <typename Input_vertex_range>
Simplex_ptr get_key(const Input_vertex_range &vertex_range);

template <typename Input_vertex_range>
std::vector<Simplex> facets(const Input_vertex_range &vertex_range);

template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1 &vertex_range1, const Input_vertex_range2 &vertex_range2);


/* sigma must not have neither face nor coface in the complex */
template <typename Input_vertex_range>
void SAL::insert_critical_simplex(const Input_vertex_range &vertex_range, Filtration_value f){
    Simplex_ptr sptr = get_key(vertex_range);
    sptr->filtration = f;
    for(const Vertex& v : vertex_range){
        if(!t0.count(v)) t0.emplace(v, Simplex_ptr_set());
        t0.at(v).emplace(sptr);
    }
}

template <typename Input_vertex_range>
void SAL::insert_simplex(const Input_vertex_range &vertex_range, Filtration_value f){
    if(membership(vertex_range)) return;
    bool all_facets_crit = true;
    for(const Simplex& facet : facets(vertex_range))
        if(!criticality(facet) || get_filtration(facet)!=f) all_facets_crit=false;
    if(all_facets_crit)
        for(const Simplex& facet : facets(vertex_range))
            erase_critical(get_key(facet));
    else
        for(const Vertex& v : vertex_range)
            //Copy constructor needed because the set is modified
            if(t0.count(v))  for(const Simplex_ptr& fptr : Simplex_ptr_set(t0.at(v)))
                if(included(*fptr,vertex_range) && get_filtration(*fptr) >= f) erase_critical(fptr); // We erase all the maximal faces of the simplex
    insert_critical_simplex(vertex_range);
}

template <typename Input_vertex_range>
void SAL::remove_simplex(const Input_vertex_range &vertex_range){
    if(vertex_range.begin()==vertex_range.end())
        t0.clear();
    else {
        const Vertex& v = best_index(vertex_range);
        //Copy constructor needed because the set is modified
        if(t0.count(v)) for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
            if(included(vertex_range, *sptr)){
                erase_critical(sptr);
                for(const Simplex& f : facets(vertex_range))
                    if(!membership(f)) insert_critical_simplex(f); // We add the facets which are new maximal simplices
            }
    }
}

template <typename Input_vertex_range>
Simplex_ptr SAL::find(const Input_vertex_range &vertex_range) const{
    if(t0.size()==0) return null_simplex();
    const Vertex& v = best_index(vertex_range);
    if(!t0.count(v))  return null_simplex();
    if(criticality(vertex_range)) return *(t0.at(v).find(get_key(vertex_range)));
    for(const Simplex_ptr& sptr : t0.at(v))
        if(included(vertex_range, *sptr)) return sptr;
    return null_simplex();
}

template <typename Input_vertex_range>
bool SAL::membership(const Input_vertex_range &vertex_range) const{
    return find(vertex_range) != null_simplex();
}

template <typename Input_vertex_range>
bool SAL::criticality(const Input_vertex_range &vertex_range) const{
    const Vertex& v =  best_index(vertex_range);
    if(!t0.count(v)) return false;
    return t0.at(v).count(get_key(vertex_range));
}

template <typename Input_vertex_range>
Simplex_ptr_set SAL::critical_cofaces(const Input_vertex_range &vertex_range) const{
    Simplex_ptr_set cofaces;
    if(criticality(vertex_range))
        cofaces.emplace(get_key(vertex_range));
    else if(vertex_range.begin()==vertex_range.end())
        for(const auto& kv : t0)
            for(const Simplex_ptr& sptr : kv.second) //kv.second is a Simplex_ptr_set
                cofaces.emplace(sptr);
    else {
        const Vertex& v = best_index(vertex_range);
        if(t0.count(v)) for(const Simplex_ptr& sptr : t0.at(v))
            if(included(vertex_range, *sptr)) cofaces.emplace(sptr);
    }
    return cofaces;
}

template <typename Input_vertex_range>
SAL::Filtration_value SAL::get_filtration(const Input_vertex_range &vertex_range) const{
    return find(vertex_range)->filtration;
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
        erase_critical(sptr);
        sigma.erase(d);
        sigma.insert(k);
        insert_simplex(sigma);
    }
    return k;
}

/* /!\ No facets added */
inline void SAL::erase_critical(const Simplex_ptr& sptr){
    if(sptr->size()==0){
        t0.at(default_vertex).erase(sptr);
        if(t0.at(default_vertex).size()==0) t0.erase(default_vertex);
    }
    else for(const Vertex& v : *sptr){
        t0.at(v).erase(sptr);
        if(t0.at(v).size()==0) t0.erase(v);
    }
}

template <typename Input_vertex_range>
Vertex SAL::best_index(const Input_vertex_range &vertex_range) const{
    std::size_t min = std::numeric_limits<size_t>::max();
    Vertex arg_min = default_vertex;
    for(const Vertex& v : vertex_range)
        if(!t0.count(v)) return v;
        else if(t0.at(v).size() < min)
            min = t0.at(v).size(), arg_min = v;
    return arg_min;
}

std::size_t SAL::num_simplices() const{
    return critical_cofaces(Simplex()).size(); //not efficient !!
}

template <typename Input_vertex_range>
bool SAL::all_facets_inside(const Input_vertex_range &vertex_range) const{
    Vertex v = best_index(vertex_range);
    if(!t0.count(v))  return false;
    Simplex sigma(vertex_range);
    Simplex f = sigma; f.erase(v);
    if(!membership(f)) return false;
    std::unordered_set<Vertex> facets_inside;
    for(const Simplex_ptr& sptr : t0.at(v))
        for(const Vertex& w : f){
            Simplex g = sigma; g.erase(w);
            if(included(g, *sptr)) facets_inside.insert(w);
        }
    return facets_inside.size() == sigma.size() - 1;
}

Simplex_ptr_set SAL::candidates() const{
    Simplex_ptr_set c;
    std::unordered_map<Simplex_ptr, std::vector<Vertex>, Sptr_hash, Sptr_equal> facets_to_max;
    for(const auto& kv : t0)
        for(const Simplex_ptr& sptr_ : kv.second){
            Simplex sigma (*sptr_);
            sigma.erase(kv.first);
            auto sptr = get_key(sigma);
            if(!facets_to_max.count(sptr)) facets_to_max.emplace(sptr, std::vector<Vertex>());
            facets_to_max.at(sptr).emplace_back(kv.first);
        }
    for(const auto& kv : facets_to_max){
        std::unordered_set<Vertex> facets(kv.second.begin(), kv.second.end());
        for(Vertex v : kv.second){
            facets.erase(v);
            for(Vertex w : facets){
                Simplex sigma(*(kv.first));
                sigma.insert(v);
                sigma.insert(w);
                if(all_facets_inside(sigma))
                    c.emplace(get_key(sigma));
            }
            facets.emplace(v);
        }
    }
    return c;
}

Simplex_ptr SAL::null_simplex() const{
    return Simplex_ptr();
}

Simplex SAL::simplex_vertex_range (const Simplex_ptr& sptr) const{
    return *sptr;
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

template <typename Input_vertex_range>
std::vector<Simplex> facets(const Input_vertex_range &vertex_range){
    std::vector<Simplex> facets;
    Simplex f(vertex_range);
    for(const Vertex& v : vertex_range){
        f.erase(v);
        facets.emplace_back(f);
        f.insert(v);
    }
    return facets;
}

template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1 &vertex_range1, const Input_vertex_range2 &vertex_range2){
    Simplex s2(vertex_range2);
    for(const Vertex& v : vertex_range1)
        if(!s2.count(v)) return false;
    return true;
}
    
template <typename Input_vertex_range>
Simplex_ptr get_key(const Input_vertex_range &vertex_range){
    Simplex s(vertex_range);
    return std::make_shared<Simplex>(s);
}



} //namespace Gudhi

#endif /* SAL_H */
