#ifndef TOWER_CONVERTER_H
#define TOWER_CONVERTER_H

#include <string>
#include <iostream>
#include <vector>

namespace Gudhi {
namespace tmp_package_name {

template<class ComplexStructure>
class Tower_converter
{
public:
    using vertex = double;
    using index = double;
    using simplex_base = std::vector<vertex>;

    enum type : int {INCLUSION, CONTRACTION, COMMENT};

    Tower_converter(std::string outputFileName);
    ~Tower_converter();

    bool add_insertion(simplex_base *simplex, double timestamp);
    bool add_contraction(vertex v, vertex u, double timestamp);

private:
    ComplexStructure complex_;
    std::ofstream *outputStream_;
};

template<class ComplexStructure>
Tower_converter::Tower_converter(std::string outputFileName)
{
    outputStream_ = new std::ofstream(outputFileName);
}

template<class ComplexStructure>
Tower_converter::~Tower_converter()
{
    delete outputStream_;
}

}
}

#endif // TOWER_CONVERTER_H
