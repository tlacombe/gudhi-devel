/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  INRIA Sophia Antipolis-Méditerranée (France)
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
//#include <list>
//#include <queue>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <cmath>
//#include <algorithm>
//#include <iomanip>
//#include <functional>

namespace Gudhi {
namespace tmp_package_name {

class Hash_complex
{
public:
    using vertex = double;
    using index = double;
    using simplex_base = std::vector<vertex>;

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

    bool insert_simplex(simplex_base *numVertices);
    bool remove_simplex(simplex_base *simplex);
    //vertex get_smallest_closed_star(vertex s1, vertex s2, std::vector<simplex_base*> *closedStar);  //returns vertex with smallest closed star; closedStar is ordered
                                                                                                    //and simplices in closedStar are independent of the ones in Complex
    //index get_boundary(simplex_base *simplex, std::vector<index> *boundary);
    //double get_size() const;

    /* Other */

    //index contract_vertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers);
    //bool exists(simplex_base *simplex);
    //void print(std::string outputFileName);
    //double get_max_size() const;

    //double get_number_of_vertices() const;
    //bool is_saturated();
    //vertex get_vertex_number(double num) const override;

private:
    std::unordered_map<double, vertex> *vertices_;
    double maxIndex_;
    double maxSize_;
    int maxDim_;
	std::unordered_map<
			std::pair<simplex_base*, int>*,
			Simplex*,
            Key_hasher,
            Simplices_equals
        > *simplices_;

    //index contract_vertex_to(std::vector<simplex_base*> *acs, simplex_base *vertexToRemoveAsSimplex, double timestamp,
    //					  std::vector<std::vector<index>*> *boundaries, std::vector<index> *&insertionNumbers);
    //bool intern_insert_simplex(simplex_base *vs, double timestamp);
    //void delete_simplex(simplex_base *simplex, std::vector<index> *&insertionNumbers);
    //int get_vertex_index(simplex_base *simplex, vertex v);
    //simplex_base* get_extended_simplex(simplex_base *simplex, vertex v);
    //vertex get_smallest_active_closed_star(simplex_base *v, simplex_base *u,
    //						   std::vector<simplex_base *> *acsActive, std::vector<simplex_base *> *acsInactive);
    //vertex get_smallest_active_star(simplex_base *v, simplex_base *u, std::queue<simplex_base*> *qv, std::queue<simplex_base*> *qu);
    //void preprocess_active_closed_star(vertex v, std::vector<simplex_base *> *&acsInactive,
    //								std::vector<simplex_base*> *&acsActive, std::vector<simplex_base*> *&acs);
    //index intern_get_boundary(simplex_base *simplex, std::vector<index> *boundary);
    //void stream_simplex(Simplex *vs, double timestamp);
    //void stream_inactivity(index insertionNumber);
};

Hash_complex::Hash_complex() : maxIndex_(-1), maxSize_(0), maxDim_(0)
{
    vertices_ = new std::unordered_map<double, vertex>();

    simplices_ = new std::unordered_map<
            std::pair<simplex_base*, int>*,
            Simplex*,
            Key_hasher,
            Simplices_equals
            >();
}

Hash_complex::~Hash_complex()
{
    for (auto it = simplices_->begin(); it != simplices_->end(); ++it){
        delete it->first;
        delete it->second;
    }
    delete simplices_;
    delete vertices_;
}

bool Hash_complex::insert_simplex(Hash_complex::simplex_base *numVertices)
{
    simplex_base *vs = new simplex_base(numVertices);
    std::pair<simplex_base*,int> *p = new std::pair<simplex_base*,int>(vs, -1);
    Simplex *splx = new Simplex(maxIndex_ + 1, vs);

    if (simplices_->emplace(p, splx).second == false) {
        p->first = NULL;
        delete splx;
        delete p;
        return false;
    }

    maxIndex_++;
    if ((int)vs->size() - 1 > maxDim_) maxDim_ = vs->size() - 1;
    if (maxSize_ < (double)simplices_->size()) maxSize_ = simplices_->size();

    if (vs->size() <= 1) return true;

    for (simplex_base::size_type i = 0; i < vs->size(); i++){
        p->second = i;
        simplices_->at(p)->add_cofacet(splx, vs->at(i));
    }
    p->second = -1;

    return true;
}

bool Hash_complex::remove_simplex(Hash_complex::simplex_base *simplex)
{
    std::pair<simplex_base*,int> p(simplex, -1);
    Simplex *splx = simplices_->at(&p);
    std::unordered_map<vertex, Simplex*> *cofacets = splx->get_cofacets();

    if (simplex->size() > 1){
        for (int i = 0; i < (int)simplex->size(); i++){
            p.second = i;
            simplices_->at(&p)->get_cofacets()->erase(simplex->at(i));
        }
    }
    p.second = -1;

    while (!cofacets->empty()) remove_simplex(cofacets->begin()->second->get_vertices());

    std::pair<simplex_base*,int> *tmp = simplices_->find(&p)->first;
    simplices_->erase(&p);
    delete tmp;
    delete splx;
}

Hash_complex::index Hash_complex::Simplex::get_insertion_num() const
{
    return insertionNum_;
}

void Hash_complex::Simplex::set_insertion_num(const index &insertionNum)
{
    insertionNum_ = insertionNum;
}

void Hash_complex::Simplex::add_cofacet(Hash_complex::Simplex *coface, Hash_complex::vertex v)
{
    cofacets_->emplace(v, coface);
}

std::unordered_map<Hash_complex::vertex, Hash_complex::Simplex*>* Hash_complex::Simplex::get_cofacets() const
{
    return cofacets_;
}

Hash_complex::simplex_base *Hash_complex::Simplex::get_vertices() const
{
    return vertices_;
}

}
}

#endif // COMPLEX_H
