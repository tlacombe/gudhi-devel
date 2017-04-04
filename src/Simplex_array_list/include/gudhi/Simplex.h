#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <unordered_set>

namespace Gudhi {

class Simplex {

public :
    typedef double Filtration_t;
    typedef std::size_t Vertex;

    Simplex();
    template <typename Input_vertex_range>
    Simplex(Input_vertex_range vertex_range);

    std::unordered_set<Vertex>::const_iterator cbegin() const;
    std::unordered_set<Vertex>::const_iterator cend() const;
    std::unordered_set<Vertex>::const_iterator begin() const;
    std::unordered_set<Vertex>::const_iterator end() const;
    std::unordered_set<Vertex>::iterator begin();
    std::unordered_set<Vertex>::iterator end();
    bool count(const Vertex v) const;
    void erase(const Vertex v);
    void insert(const Vertex v);
    std::size_t size() const;

    Filtration_t filtration;

private:
    std::unordered_set<Vertex> vertices;
};

Simplex::Simplex()
    : vertices()
{}

template <typename Input_vertex_range>
Simplex::Simplex(Input_vertex_range vertex_range)
    : vertices(vertex_range.begin(),vertex_range.end())
{}

std::unordered_set<Simplex::Vertex>::const_iterator Simplex::cbegin() const{
    return vertices.cbegin();
}

std::unordered_set<Simplex::Vertex>::const_iterator Simplex::cend() const{
    return vertices.cend();
}

std::unordered_set<Simplex::Vertex>::const_iterator Simplex::begin() const{
    return vertices.begin();
}

std::unordered_set<Simplex::Vertex>::const_iterator Simplex::end() const{
    return vertices.end();
}

  
std::unordered_set<Simplex::Vertex>::iterator Simplex::begin(){
    return vertices.begin();
}

std::unordered_set<Simplex::Vertex>::iterator Simplex::end(){
    return vertices.end();
}
  
bool Simplex::count(const Vertex v) const{
    return vertices.count(v);
}

void Simplex::erase(const Vertex v){
    vertices.erase(v);
}

void Simplex::insert(const Vertex v){
    vertices.insert(v);
}

std::size_t Simplex::size() const{
    return vertices.size();
}

} //namespace Gudhi

#endif /* SIMPLEX_H */
