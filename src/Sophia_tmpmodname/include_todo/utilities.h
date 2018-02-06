#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <string>

#include "tower_converter.h"
#include "Persistent_homology/persistence.h"

namespace Gudhi {
namespace tmp_package_name {

template<class ComplexStructure>
typename Tower_converter<ComplexStructure>::operationType read_operation(std::string *line, std::vector<double> *vertices, double *timestamp)
{
    using TC = Tower_converter<ComplexStructure>;
    typename TC::operationType type;
    vertices->clear();
    double num;

    size_t next = line->find_first_not_of(' ', 0);
    size_t current = next;
    next = line->find_first_of(' ', current);
    if (next == std::string::npos) return TC::COMMENT;
    if (line->substr(current, next - current) == "i") type = TC::INCLUSION;
    else if (line->substr(current, next - current) == "c") type = TC::CONTRACTION;
    else if (line->substr(current, next - current) == "#") return TC::COMMENT;
    else {
        *timestamp = stod(line->substr(current, next - current));
        next = line->find_first_not_of(' ', next + 1);
        current = next;
        next = line->find_first_of(' ', current);
        if (next == std::string::npos) {
            std::cout << "Operation syntaxe error in file.\n";
            exit(0);
        }
        if (line->substr(current, next - current) == "i") type = TC::INCLUSION;
        else if (line->substr(current, next - current) == "c") type = TC::CONTRACTION;
        else if (line->substr(current, next - current) == "#") return TC::COMMENT;
        else {
            std::cout << "Operation syntaxe error in file.\n";
            exit(0);
        }
    }

    next = line->find_first_not_of(' ', next + 1);
    while (next != std::string::npos){
        current = next;
        next = line->find_first_of(' ', current);
        num = stod(line->substr(current, next - current));
        vertices->push_back(num);
        if (next != std::string::npos) next = line->find_first_not_of(' ', next + 1);
    }

    return type;
}

template<class ComplexStructure>
std::ifstream& operator>>(std::ifstream& file, Tower_converter<ComplexStructure>& tc)
{
    using TC = Tower_converter<ComplexStructure>;
    std::string line;

    if (file.is_open()){
        std::vector<double> vertices;
        double timestamp = -1;
        double defaultTimestamp = 0;
        while (getline(file, line, '\n')){
            typename TC::operationType type = read_operation<ComplexStructure>(&line, &vertices, &timestamp);
            if (timestamp != -1) defaultTimestamp = timestamp;

            if (type == TC::INCLUSION){
                if (tc.add_insertion(&vertices, defaultTimestamp)) defaultTimestamp++;
            } else if (type == TC::CONTRACTION) {
                tc.add_contraction(vertices.at(0), vertices.at(1), defaultTimestamp);
                defaultTimestamp++;
            }

            timestamp = -1;
        }
        file.close();
    } else {
        std::cout << "Unable to open input file.\n";
        file.setstate(std::ios::failbit);
    }

    return file;
}

template<class ComplexStructure>
std::ifstream& operator>>(std::ifstream& file, Persistence<ComplexStructure>& pers)
{
    using TC = Tower_converter<ComplexStructure>;
    std::string line;

    if (file.is_open()){
        std::vector<double> vertices;
        double timestamp = -1;
        double defaultTimestamp = 0;
        while (getline(file, line, '\n')){
            typename TC::operationType type = read_operation<ComplexStructure>(&line, &vertices, &timestamp);
            if (timestamp != -1) defaultTimestamp = timestamp;

            if (type == TC::INCLUSION){
                if (pers.add_insertion(&vertices, defaultTimestamp)) defaultTimestamp++;
            } else if (type == TC::CONTRACTION) {
                pers.add_contraction(vertices.at(0), vertices.at(1), defaultTimestamp);
                defaultTimestamp++;
            }

            timestamp = -1;
        }
        pers.finalize_reduction();
        file.close();
    } else {
        std::cout << "Unable to open input file.\n";
        file.setstate(std::ios::failbit);
    }

    return file;
}

}
}

#endif // UTILITIES_H
