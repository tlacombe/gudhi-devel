#ifndef STW_H
#define STW_H

#include <gudhi/SAL.h>
#include <gudhi/Simplex_tree.h>

namespace Gudhi {

typedef Simplex typeVectorVertex;
typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

class STW {
    
public:
    void add(const Simplex& tau);
    bool membership(const Simplex& tau);
    Vertex contraction(const Vertex x, const Vertex y);
    std::size_t size();

private:
    Simplex_tree<> simplexTree;
    void erase_max(const Simplex& sigma);
};

void STW::add(const Simplex& tau){
    simplexTree.insert_simplex_and_subfaces(tau);
}

bool STW::membership(const Simplex& tau) {
    return simplexTree.find(tau) != simplexTree.null_simplex();
}

void STW::erase_max(const Simplex& sigma){
    if(membership(sigma))
        simplexTree.remove_maximal_simplex(simplexTree.find(sigma));
}

Vertex STW::contraction(const Vertex x, const Vertex y){
    Simplex sx; sx.insert(x);
    auto hx = simplexTree.find(sx);
    if(hx != simplexTree.null_simplex())
        for(auto h : simplexTree.cofaces_simplex_range(hx,0)){
            auto sr = simplexTree.simplex_vertex_range(h);
            Simplex sigma(sr.begin(),sr.end());
            erase_max(sigma);
            sigma.erase(x);
            sigma.insert(y);
            add(sigma);
        }
    return y;
}

std::size_t STW::size(){
    return simplexTree.num_simplices();
}


} //namespace Gudhi

#endif /* STW_H */
