/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  TU Graz (Austria)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMPLEX_H
#define COMPLEX_H

#include <unordered_map>
#include <vector>
#include <queue>
#include <algorithm>

namespace Gudhi {
namespace tower_to_filtration {

class Hash_complex
{
public:
    using vertex = long long;
    using index = long long;
    using simplex_base = std::vector<vertex>;
    using size_type = long long;

    Hash_complex();
    ~Hash_complex();

	class Simplex
	{
	public:
        Simplex(index num, simplex_base *vertices);
		~Simplex();

        index get_insertion_num() const;
        void set_insertion_num(const index &insertionNum);

        void add_cofacet(Simplex *coface, vertex v);
        std::unordered_map<vertex, Simplex*>* get_cofacets() const;

        simplex_base *get_vertices() const;

    private:
        index insertionNum_;
        std::unordered_map<vertex, Simplex*> *cofacets_;
        simplex_base* vertices_;
	};

    struct Key_hasher {
		std::size_t operator()(const std::pair<simplex_base*, int> *k) const
		{
			std::size_t seed;
			if (k->second < 0) seed = k->first->size();
			else seed = k->first->size() - 1;

			for (int i = 0; i < (int)k->first->size(); i++) {
				if (i != k->second) seed ^= (std::size_t)(k->first->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};

    struct Simplices_equals : std::binary_function<const std::pair<simplex_base*, int>*, const std::pair<simplex_base*, int>*, bool> {
		bool operator()(const std::pair<simplex_base*, int> *s1, const std::pair<simplex_base*, int> *s2) const
		{
			const std::pair<simplex_base*, int> *key;
			const std::pair<simplex_base*, int> *inMap;
			simplex_base::size_type size;

			if (s1->second > -1) {
				key = s1;
				inMap = s2;
			} else {
				key = s2;
				inMap = s1;
			}

			if (key->second < 0) size = key->first->size();
			else size = key->first->size() - 1;
			if (size != inMap->first->size()) return false;
			int j = 0;
			for (simplex_base::size_type i = 0; i < size; i++){
				if (j == key->second) j++;
				if (key->first->at(j) != inMap->first->at(i)) {
					return false;
				}
				j++;
			}
			return true;
		}
	};

    /* Important for tower_converter.h */

    bool insert_simplex(simplex_base &numVertices);
    bool remove_simplex(simplex_base &simplex);
    bool remove_simplex(simplex_base &simplex, std::vector<index> *removedIndices);
    vertex get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_base *> *closedStar);   //returns vertex with smallest closed star; closedStar is ordered
                                                                                                    //and simplices in closedStar are independent of the ones in Complex
    index get_boundary(simplex_base &simplex, std::vector<index> *boundary);    //ordered by increasing insertion numbers
    size_type get_size() const;
    index get_max_index() const;

    /* Other */

    //index contract_vertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers);
    size_type get_max_size() const;
    int get_max_dimension() const;
    void get_cofaces(simplex_base &simplex, std::vector<simplex_base*> *cofaces);
    void get_cofacets(simplex_base &simplex, std::vector<simplex_base *> *cofacets);
    bool contains(simplex_base &simplex);

private:
    index maxIndex_;
    size_type maxSize_;
    int maxDim_;
	std::unordered_map<
			std::pair<simplex_base*, int>*,
			Simplex*,
            Key_hasher,
            Simplices_equals
        > *simplices_;

    vertex get_smallest_star(vertex v, vertex u, std::queue<simplex_base*> *qv, std::queue<simplex_base*> *qu);
    int get_vertex_index(simplex_base *simplex, vertex v);
    //index contract_vertex_to(std::vector<simplex_base*> *acs, simplex_base *vertexToRemoveAsSimplex, double timestamp,
    //					  std::vector<std::vector<index>*> *boundaries, std::vector<index> *&insertionNumbers);
    //simplex_base* get_extended_simplex(simplex_base *simplex, vertex v);
    //void preprocess_active_closed_star(vertex v, std::vector<simplex_base *> *&acsInactive,
    //								std::vector<simplex_base*> *&acsActive, std::vector<simplex_base*> *&acs);
};

inline Hash_complex::Hash_complex() : maxIndex_(-1), maxSize_(0), maxDim_(0)
{
    simplices_ = new std::unordered_map<
            std::pair<simplex_base*, int>*,
            Simplex*,
            Key_hasher,
            Simplices_equals
            >();
}

inline Hash_complex::~Hash_complex()
{
    for (auto it = simplices_->begin(); it != simplices_->end(); ++it){
        delete it->first;
        delete it->second;
    }
    delete simplices_;
}

inline bool Hash_complex::insert_simplex(simplex_base &numVertices)
{
    simplex_base *vs = new simplex_base(numVertices);
    std::pair<simplex_base*,int> *p = new std::pair<simplex_base*,int>(vs, -1);
    Simplex *splx = new Simplex(maxIndex_ + 1, vs);

    if (simplices_->emplace(p, splx).second == false) {
	p->first = nullptr;
        delete splx;
        delete p;
        return false;
    }

    maxIndex_++;
    if ((int)vs->size() - 1 > maxDim_) maxDim_ = vs->size() - 1;
    if (maxSize_ < (size_type)simplices_->size()) maxSize_ = simplices_->size();

    if (vs->size() <= 1) return true;

    for (simplex_base::size_type i = 0; i < vs->size(); i++){
        p->second = i;
        simplices_->at(p)->add_cofacet(splx, vs->at(i));
    }
    p->second = -1;

    return true;
}

inline bool Hash_complex::remove_simplex(simplex_base &simplex)
{
    return remove_simplex(simplex, nullptr);
}

inline bool Hash_complex::remove_simplex(simplex_base &simplex, std::vector<index> *removedIndices)
{
    std::pair<simplex_base*,int> p(&simplex, -1);
    if (simplices_->find(&p) == simplices_->end()) return false;
    Simplex *splx = simplices_->at(&p);
    std::unordered_map<vertex, Simplex*> *cofacets = splx->get_cofacets();

    if (simplex.size() > 1){
	for (int i = 0; i < (int)simplex.size(); i++){
            p.second = i;
	    simplices_->at(&p)->get_cofacets()->erase(simplex.at(i));
        }
    }
    p.second = -1;

    while (!cofacets->empty()) remove_simplex(*(cofacets->begin()->second->get_vertices()), removedIndices);

    std::pair<simplex_base*,int> *tmp = simplices_->find(&p)->first;
    if (removedIndices != nullptr) removedIndices->push_back(splx->get_insertion_num());
    simplices_->erase(&p);
    delete tmp;
    delete splx;
    return true;
}

inline Hash_complex::vertex Hash_complex::get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_base *> *closedStar)
{
    std::queue<simplex_base*> qv;
    std::queue<simplex_base*> qu;
    simplex_base *s;
    std::pair<simplex_base*,int> p;

    if (get_smallest_star(v, u, &qv, &qu) == v){
        while (!qv.empty()){
            s = qv.front();
            qv.pop();
            p.first = s;
            p.second = get_vertex_index(s, v);
            if (s->size() > 1){
                closedStar->push_back(new simplex_base(*(simplices_->at(&p)->get_vertices())));
            }
            p.second = -1;
            closedStar->push_back(new simplex_base(*s));
        }

        return v;
    } else {
        while (!qu.empty()){
            s = qu.front();
            qu.pop();
            p.first = s;
            p.second = get_vertex_index(s, u);
            if (s->size() > 1){
                closedStar->push_back(new simplex_base(*(simplices_->at(&p)->get_vertices())));
            }
            p.second = -1;
            closedStar->push_back(new simplex_base(*s));
        }

        return u;
    }
}

inline Hash_complex::index Hash_complex::get_boundary(simplex_base &simplex, std::vector<index> *boundary)
{
    std::pair<simplex_base*,int> p(&simplex, -1);

    if (simplex.size() == 1) return simplices_->at(&p)->get_insertion_num();

    for (simplex_base::size_type i = 0; i < simplex.size(); i++){
        p.second = i;
        boundary->push_back(simplices_->at(&p)->get_insertion_num());
    }
    std::sort(boundary->begin(), boundary->end());

    p.second = -1;
    return simplices_->at(&p)->get_insertion_num();
}

inline Hash_complex::size_type Hash_complex::get_size() const
{
    return simplices_->size();
}

inline Hash_complex::index Hash_complex::get_max_index() const
{
    return maxIndex_;
}

inline Hash_complex::size_type Hash_complex::get_max_size() const
{
    return maxSize_;
}

inline int Hash_complex::get_max_dimension() const
{
    return maxDim_;
}

inline void Hash_complex::get_cofaces(simplex_base &simplex, std::vector<simplex_base *> *cofaces)
{
    auto hash = [](simplex_base* const& k) {
        std::size_t seed = k->size();
        for (simplex_base::size_type i = 0; i < k->size(); i++) {
            seed ^= (std::size_t)(k->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    };
    auto comp = [](const simplex_base* s1, const simplex_base* s2) {
        simplex_base::size_type size = s1->size();
        if (size != s2->size()) return false;
        for (simplex_base::size_type i = 0; i < size; i++){
            if (s1->at(i) != s2->at(i)) return false;
        }
        return true;
    };

    std::unordered_map<simplex_base*, bool, decltype(hash), decltype(comp)> visited(100, hash, comp);
    std::pair<simplex_base*,int> p(&simplex, -1);
    std::unordered_map<vertex, Simplex*> *cofacets;
    std::unordered_map<vertex, Simplex*>::iterator it;
    std::queue<simplex_base*> q;

    q.push(&simplex);

    while (!q.empty()){
        p.first = q.front();
        q.pop();
        cofacets = simplices_->at(&p)->get_cofacets();
        for (it = cofacets->begin(); it != cofacets->end(); it++){
            p.first = it->second->get_vertices();
            if (visited.find(p.first) == visited.end()){
                visited.emplace(p.first, true);
                cofaces->push_back(new simplex_base(*(p.first)));
                q.push(p.first);
            }
        }
    }
}

inline void Hash_complex::get_cofacets(simplex_base &simplex, std::vector<simplex_base *> *cofacets)
{
    std::pair<simplex_base*,int> p(&simplex, -1);
    std::unordered_map<vertex, Simplex*> *mapedCofacets = simplices_->at(&p)->get_cofacets();
    for (auto it = mapedCofacets->begin(); it != mapedCofacets->end(); it++){
        cofacets->push_back(new simplex_base(*(it->second->get_vertices())));
    }
}

inline bool Hash_complex::contains(simplex_base &simplex)
{
    std::pair<simplex_base*,int> p(&simplex, -1);
    return (simplices_->find(&p) != simplices_->end());
}

inline Hash_complex::vertex Hash_complex::get_smallest_star(vertex v, vertex u, std::queue<simplex_base *> *qv, std::queue<simplex_base *> *qu)
{
    auto hash = [](simplex_base* const& k) {
        std::size_t seed = k->size();
        for (simplex_base::size_type i = 0; i < k->size(); i++) {
            seed ^= (std::size_t)(k->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    };
    auto comp = [](const simplex_base* s1, const simplex_base* s2) {
        simplex_base::size_type size = s1->size();
        if (size != s2->size()) return false;
        for (simplex_base::size_type i = 0; i < size; i++){
            if (s1->at(i) != s2->at(i)) return false;
        }
        return true;
    };

    std::unordered_map<simplex_base*, bool, decltype(hash), decltype(comp)> visitedV(100, hash, comp);
    std::unordered_map<simplex_base*, bool, decltype(hash), decltype(comp)> visitedU(100, hash, comp);
    std::pair<simplex_base*,int> p;
    std::unordered_map<vertex, Simplex*> *cofacetsV;
    std::unordered_map<vertex, Simplex*> *cofacetsU;
    std::unordered_map<vertex, Simplex*>::iterator vit;
    std::unordered_map<vertex, Simplex*>::iterator uit;
    std::queue<simplex_base*> tv;
    std::queue<simplex_base*> tu;
    simplex_base sv(1, v);
    simplex_base su(1, u);

    p.first = &sv;
    p.second = -1;
    qv->push(simplices_->at(&p)->get_vertices());
    cofacetsV = simplices_->at(&p)->get_cofacets();
    vit = cofacetsV->begin();

    p.first = &su;
    qu->push(simplices_->at(&p)->get_vertices());
    cofacetsU = simplices_->at(&p)->get_cofacets();
    uit = cofacetsU->begin();

    if (vit != cofacetsV->end() && vit->first == u) vit++;
    if (uit != cofacetsU->end() && uit->first == v) uit++;

    while (vit != cofacetsV->end() && uit != cofacetsU->end()) {
        p.first = vit->second->get_vertices();
        visitedV.emplace(p.first, true);
        qv->push(p.first);
        tv.push(p.first);

        p.first = uit->second->get_vertices();
        visitedU.emplace(p.first, true);
        qu->push(p.first);
        tu.push(p.first);

        vit++;
        while ((vit == cofacetsV->end() && !tv.empty()) || (vit != cofacetsV->end() && (vit->first == u || visitedV.find(vit->second->get_vertices()) != visitedV.end()))){
            if (vit == cofacetsV->end() && !tv.empty()){
                p.first = tv.front();
                tv.pop();
                cofacetsV = simplices_->at(&p)->get_cofacets();
                vit = cofacetsV->begin();
            }
            if (vit != cofacetsV->end() && (vit->first == u || visitedV.find(vit->second->get_vertices()) != visitedV.end())) vit++;
        }

        uit++;
        while ((uit == cofacetsU->end() && !tu.empty()) || (uit != cofacetsU->end() && (uit->first == v || visitedU.find(uit->second->get_vertices()) != visitedU.end()))){
            if (uit == cofacetsU->end() && !tu.empty()){
                p.first = tu.front();
                tu.pop();
                cofacetsU = simplices_->at(&p)->get_cofacets();
                uit = cofacetsU->begin();
            }
            if (uit != cofacetsU->end() && (uit->first == v || visitedU.find(uit->second->get_vertices()) != visitedU.end())) uit++;
        }
    }

    if (uit == cofacetsU->end()) return u;
    else return v;
}

inline int Hash_complex::get_vertex_index(simplex_base *simplex, vertex v)
{
    int i = 0;
    while (i < (int)simplex->size() && simplex->at(i) != v){
        i++;
    }
    if (i == (int)simplex->size()) return -1;
    return i;
}

inline Hash_complex::Simplex::Simplex(index num, simplex_base *vertices) : insertionNum_(num), vertices_(vertices)
{
    cofacets_ = new std::unordered_map<vertex, Simplex*>();
}

inline Hash_complex::Simplex::~Simplex()
{
    delete cofacets_;
    delete vertices_;
}

inline Hash_complex::index Hash_complex::Simplex::get_insertion_num() const
{
    return insertionNum_;
}

inline void Hash_complex::Simplex::set_insertion_num(const index &insertionNum)
{
    insertionNum_ = insertionNum;
}

inline void Hash_complex::Simplex::add_cofacet(Simplex *coface, vertex v)
{
    cofacets_->emplace(v, coface);
}

inline std::unordered_map<Hash_complex::vertex, Hash_complex::Simplex*>* Hash_complex::Simplex::get_cofacets() const
{
    return cofacets_;
}

inline Hash_complex::simplex_base *Hash_complex::Simplex::get_vertices() const
{
    return vertices_;
}

}
}

#endif // COMPLEX_H
