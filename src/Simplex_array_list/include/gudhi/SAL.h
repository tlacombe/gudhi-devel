#ifndef SAL_H
#define SAL_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <limits>

namespace Gudhi {

typedef std::size_t Vertex;
typedef std::unordered_set<Vertex> Simplex;
typedef std::shared_ptr<Simplex> Simplex_ptr;
struct Sptr_hash{ std::size_t operator()(const Simplex_ptr& s) const; };
struct Sptr_equal{ std::size_t operator()(const Simplex_ptr& a, const Simplex_ptr& b) const; };
typedef std::unordered_set<Simplex_ptr, Sptr_hash, Sptr_equal> Simplex_ptr_set;

Vertex default_vertex = 0;
Simplex_ptr null_simplex_ptr = Simplex_ptr();

template <typename Input_vertex_range>
Simplex_ptr get_key(const Input_vertex_range &vertex_range);

template <typename Input_vertex_range>
std::vector<Simplex> facets(const Input_vertex_range &vertex_range);

template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1 &vertex_range1, const Input_vertex_range2 &vertex_range2);

class SAL {
    
public:
    template <typename Input_vertex_range>
    void insert_max(const Input_vertex_range &vertex_range);

    template <typename Input_vertex_range>
    void add(const Input_vertex_range &vertex_range);
    template <typename Input_vertex_range>
    void remove(const Input_vertex_range &vertex_range);
    
    template <typename Input_vertex_range>
    Simplex_ptr find(const Input_vertex_range &vertex_range) const;
    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range) const;
    template <typename Input_vertex_range>
    bool all_facets_inside(const Input_vertex_range &vertex_range) const;
    template <typename Input_vertex_range>
    bool maximality(const Input_vertex_range &vertex_range) const;
    template <typename Input_vertex_range>
    Simplex_ptr_set max_cofaces(const Input_vertex_range &vertex_range) const;


    Vertex contraction(const Vertex x, const Vertex y);

    std::size_t size() const;
    std::size_t num_vertices() const;
    

  // protected:
    void erase_max(const Simplex_ptr& sptr);
    template <typename Input_vertex_range>
    Vertex best_index(const Input_vertex_range &vertex_range) const;
    
    std::unordered_map<Vertex, Simplex_ptr_set> t0;
    bool max_empty_face; // Is the empty simplex a maximal face ?

};


/* sigma must not have neither face nor coface in the complex */
template <typename Input_vertex_range>
void SAL::insert_max(const Input_vertex_range &vertex_range){
    Simplex_ptr sptr = get_key(vertex_range);
    max_empty_face = (sptr->size()==0);
    for(const Vertex& v : vertex_range){
        if(!t0.count(v)) t0.emplace(v, Simplex_ptr_set());
        t0.at(v).emplace(sptr);
    }
}

template <typename Input_vertex_range>
void SAL::add(const Input_vertex_range &vertex_range){
    if(membership(vertex_range)) return;
    bool all_facets_max = true;
    for(const Simplex& f : facets(vertex_range))
        if(!maximality(f)) all_facets_max=false;
    if(all_facets_max)
        for(const Simplex& f : facets(vertex_range))
            erase_max(get_key(f));
    else
        for(const Vertex& v : vertex_range)
            //Copy constructor needed because the set is modified
            if(t0.count(v))  for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
                if(included(*sptr,vertex_range)) erase_max(sptr); // We erase all the maximal faces of the simplex
    insert_max(vertex_range);
}

template <typename Input_vertex_range>
void SAL::remove(const Input_vertex_range &vertex_range){
    if(vertex_range.begin()==vertex_range.end()){
        t0.clear();
        max_empty_face = false;
    }
    else {
        const Vertex& v = best_index(vertex_range);
        //Copy constructor needed because the set is modified
        if(t0.count(v)) for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
            if(included(vertex_range, *sptr)){
                erase_max(sptr);
                for(const Simplex& f : facets(vertex_range))
                    if(!membership(f)) insert_max(f); // We add the facets which are new maximal simplices
            }
    }
}

template <typename Input_vertex_range>
Simplex_ptr SAL::find(const Input_vertex_range &vertex_range) const{
    if(t0.size()==0) return null_simplex_ptr;
    const Vertex& v = best_index(vertex_range);
    if(!t0.count(v))  return null_simplex_ptr;
    if(maximality(vertex_range)) return *(t0.at(v).find(get_key(vertex_range)));
    for(const Simplex_ptr& sptr : t0.at(v))
        if(included(vertex_range, *sptr)) return sptr;
    return null_simplex_ptr;
}

template <typename Input_vertex_range>
bool SAL::membership(const Input_vertex_range &vertex_range) const{
    return find(vertex_range) != null_simplex_ptr;
}

template <typename Input_vertex_range>
bool SAL::all_facets_inside(const Input_vertex_range &vertex_range) const{
    Vertex v = best_index(vertex_range);
    if(!t0.count(v))  return false;
    Simplex sigma(vertex_range.begin(),vertex_range.end());
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


template <typename Input_vertex_range>
bool SAL::maximality(const Input_vertex_range &vertex_range) const{
    if(t0.size()==0 && !max_empty_face) return false; //empty complex
    const Vertex& v =  best_index(vertex_range);
    if(!t0.count(v)) return false;
    return t0.at(v).count(get_key(vertex_range));
}

template <typename Input_vertex_range>
Simplex_ptr_set SAL::max_cofaces(const Input_vertex_range &vertex_range) const{
    Simplex_ptr_set cofaces;
    if(maximality(vertex_range))
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
        erase_max(sptr);
        sigma.erase(d);
        sigma.insert(k);
        add(sigma);
    }
    return k;
}

/* /!\ No facets added */
inline void SAL::erase_max(const Simplex_ptr& sptr){
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

std::size_t SAL::size() const{
    return max_cofaces(Simplex()).size(); //not efficient
}

/* Return the number of vertices
 */
std::size_t SAL::num_vertices() const {
  return t0.size();
}

std::size_t Sptr_equal::operator()(const Simplex_ptr& s1, const Simplex_ptr& s2) const {
    if (s1->size() != s2->size()) return false;
    return included(*s1,*s2); //tests equality for same size simplices
}


std::size_t Sptr_hash::operator()(const Simplex_ptr& s) const {
    std::hash<double> h_f; //double hash is better than int hash
    size_t h = 0;
    for(const Vertex& v : *s)
        h += h_f(static_cast<double>(v));
    return h;
}

template <typename Input_vertex_range>
std::vector<Simplex> facets(const Input_vertex_range &vertex_range){
    std::vector<Simplex> facets;
    Simplex f(vertex_range.begin(),vertex_range.end());
    for(const Vertex& v : vertex_range){
        f.erase(v);
        facets.emplace_back(f);
        f.insert(v);
    }
    return facets;
}

template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1 &vertex_range1, const Input_vertex_range2 &vertex_range2){
    Simplex s2(vertex_range2.begin(),vertex_range2.end());
    for(const Vertex& v : vertex_range1)
        if(!s2.count(v)) return false;
    return true;
}
    
template <typename Input_vertex_range>
Simplex_ptr get_key(const Input_vertex_range &vertex_range){
    Simplex s(vertex_range.begin(),vertex_range.end());
    return std::make_shared<Simplex>(s);
}

} //namespace Gudhi

#endif /* SAL_H */
