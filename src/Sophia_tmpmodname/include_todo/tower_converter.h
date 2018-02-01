#ifndef TOWER_CONVERTER_H
#define TOWER_CONVERTER_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <algorithm>

namespace Gudhi {
namespace tmp_package_name {

template<class ComplexStructure>
class Tower_converter
{
public:
    using vertex = double;
    using index = double;
    using simplex_base = std::vector<vertex>;

    enum operationType : int {INCLUSION, CONTRACTION, COMMENT};
    enum streamingType : int {FACES, VERTICES};

    Tower_converter(ComplexStructure *complex, std::string outputFileName);
    ~Tower_converter();

    bool add_insertion(std::vector<double> *simplex, double timestamp);
    void add_contraction(double v, double u, double timestamp);

    double get_filtration_size() const;
    double get_tower_width() const;

private:
    ComplexStructure *complex_;
    std::unordered_map<double, vertex> *vertices_;
    std::ofstream *outputStream_;
    streamingType streamingType_;
    double filtrationSize_;
    double towerWidth_;

    void get_union(vertex v, std::vector<simplex_base*> *simplices);
    void stream_simplex(simplex_base *simplex, double timestamp);
};

template<class ComplexStructure>
Tower_converter<ComplexStructure>::Tower_converter(ComplexStructure *complex, std::string outputFileName)
    : complex_(complex), streamingType_(VERTICES), filtrationSize_(0), towerWidth_(0)
{
    outputStream_ = new std::ofstream(outputFileName);
    vertices_ = new std::unordered_map<double, vertex>();
    //streamingType_ = FACES;
}

template<class ComplexStructure>
Tower_converter<ComplexStructure>::~Tower_converter()
{
    delete outputStream_;
    delete vertices_;
}

template<class ComplexStructure>
bool Tower_converter<ComplexStructure>::add_insertion(std::vector<double> *simplex, double timestamp)
{
    simplex_base transSimplex;

    if (simplex->size() == 1){
        vertices_->emplace(simplex->at(0), simplex->at(0));
        transSimplex.push_back(simplex->at(0));
    } else {
        for (simplex_base::size_type i = 0; i < simplex->size(); i++){
            transSimplex.push_back(vertices_->at(simplex->at(i)));
        }
        std::sort(transSimplex.begin(), transSimplex.end());
    }

    if (complex_->insert_simplex(&transSimplex)) {
        stream_simplex(&transSimplex, timestamp);
        if (complex_->get_size() > towerWidth_) towerWidth_ = complex_->get_size();
        return true;
    }
    return false;
}

template<class ComplexStructure>
void Tower_converter<ComplexStructure>::add_contraction(double v, double u, double timestamp)
{
    std::vector<simplex_base*> closedStar;
    vertex tv = vertices_->at(v), tu = vertices_->at(u);
    vertex dis = complex_->get_smallest_closed_star(tv, tu, &closedStar);
    simplex_base vdis(1, dis);

    vertices_->erase(v);
    if (dis == tu){
        vertices_->at(u) = tv;
        get_union(tv, &closedStar);
    } else {
        get_union(tu, &closedStar);
    }

    for (auto it = closedStar.begin(); it != closedStar.end(); it++){
        if (complex_->insert_simplex(*it)) stream_simplex(*it, timestamp);
        delete *it;
    }
    complex_->remove_simplex(&vdis);

    if (complex_->get_size() > towerWidth_) towerWidth_ = complex_->get_size();
}

template<class ComplexStructure>
double Tower_converter<ComplexStructure>::get_filtration_size() const
{
    return filtrationSize_;
}

template<class ComplexStructure>
double Tower_converter<ComplexStructure>::get_tower_width() const
{
    return towerWidth_;
}

template<class ComplexStructure>
void Tower_converter<ComplexStructure>::get_union(vertex v, std::vector<simplex_base*> *simplices)
{
    for (auto itSimplices = simplices->begin(); itSimplices != simplices->end(); itSimplices++){
        auto itVertices = (*itSimplices)->begin();
        while (itVertices != (*itSimplices)->end() && *itVertices < v) itVertices++;
        if ((itVertices != (*itSimplices)->end() && *itVertices != v) || itVertices == (*itSimplices)->end()){
            (*itSimplices)->insert(itVertices, v);
        }
    }
}

template<class ComplexStructure>
void Tower_converter<ComplexStructure>::stream_simplex(simplex_base *simplex, double timestamp)
{
    filtrationSize_++;
    simplex_base::size_type size = simplex->size();
    *outputStream_ << std::setprecision(std::numeric_limits<double>::digits10 + 1) << (size - 1) << " ";
    if (streamingType_ == FACES){
        if (size > 1){
            std::vector<index> boundary;
            complex_->get_boundary(simplex, &boundary);
            for (std::vector<index>::size_type i = 0; i < size; i++){
                *outputStream_ << boundary.at(i) << " ";
            }
        }
    } else {
        for (simplex_base::size_type i = 0; i < size; i++){
            *outputStream_ << simplex->at(i) << " ";
        }
    }
    *outputStream_ << timestamp << "\n";
}

}
}

#endif // TOWER_CONVERTER_H
