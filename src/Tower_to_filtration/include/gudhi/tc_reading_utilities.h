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

#ifndef UTILITIES_H
#define UTILITIES_H

/** @file tc_reading_utilities.h
 * @brief Contains file reading utilities.
 */

#include <iostream>
#include <string>

#include "gudhi/tower_converter.h"
#include "gudhi/Persistent_homology/persistence.h"

namespace Gudhi {
namespace tower_to_filtration {

template<class ComplexStructure>
/**
 * @brief Reads the tower operation stored in @p line and store the corresponding vertices in @p vertices. Returns the operation type.
 *
 * For the right line format, see @ref sophiafileformat.
 *
 * @param line line te be read.
 * @param vertices pointer to an empty vector of doubles. If the operation is an inclusion, it will store the vertices of the inserted simplex in increasing order ;
 * if the operation is a contraction, it will store the two contracted vertices.
 * @param timestamp time value associated to the operation.
 * @return The operation type: #INCLUSION if it is an inclusion, #CONTRACTION if it is a contraction, or #COMMENT if it is not an operation.
 */
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
/**
 * @brief Reads @p file containing a tower and feed it to @tc to construct the corresponding filtration.
 * @param file file to be read. For the right file format, see @ref sophiafileformat.
 * @param tc instance of @ref Gudhi::tower_to_filtration::Tower_converter<ComplexStructure>.
 * @return @p file.
 */
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

template<class ComplexStructure, class ColumnType>
/**
 * @brief Reads @p file containing a tower and feed it to @p pers to compute its persistence barcode.
 * @param file file to be read. For the right file format, see @ref sophiafileformat.
 * @param pers instance of @ref Gudhi::tower_to_filtration::Persistence<ComplexStructure,ColumnType>.
 * @return @p file.
 */
std::ifstream& operator>>(std::ifstream& file, Persistence<ComplexStructure,ColumnType>& pers)
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
