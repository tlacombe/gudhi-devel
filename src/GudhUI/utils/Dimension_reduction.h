/*
 * Dimension_reduction.h
 *  Created on: Feb 6, 2015
 * This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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
 * 
 */


#ifndef DIMENSION_REDUCTION_H_
#define DIMENSION_REDUCTION_H_

#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>


template<typename Complex>
class Dimension_reduction{
	typedef typename Complex::GT GT;
	typedef typename GT::Point Point;
	typedef typename GT::Point_3 Point_3;

public:
	Dimension_reduction(const Complex& complex,int reduced_dimension=3):
		complex_(complex),reduced_dimension_(reduced_dimension){

		if(!complex.empty()){
			save_points_to_arff_file();
			run_waffle_dim_reduction();
			read_points_from_arff_file();
		}

	}


	std::vector<Point_3> reduced_points() const{
		return reduced_points_;
	}

private:

	void run_waffle_dim_reduction(){
		std::string cmd(
						"waffles_dimred isomap "+name_arff_file+
						" kdtree 14 "+std::to_string(reduced_dimension_)+" > tmp_red.arff");
//		std::string cmd(
//				"waffles_dimred breadthfirstunfolding "+name_arff_file+
//				" kdtree 14 "+std::to_string(reduced_dimension_)+" -reps 20 > tmp_red.arff");
		std::cout<<cmd<<std::endl;
		system(cmd.c_str());
	}

	void save_points_to_arff_file() const{
		std::ofstream file(name_arff_file);
		if(file.is_open()){
			file<<"@RELATION Untitled\n\n";

			auto first_point = complex_.point(*complex_.vertex_range().begin());
			int ambient_dimension = first_point.dimension();

			for(int i = 0 ; i < ambient_dimension; ++i)
				file<<"@ATTRIBUTE attr_"<<i<<"	real\n";

			file << "\n@DATA\n";
			for(auto v : complex_.vertex_range()){
				bool first = true;
				for(auto x_it = complex_.point(v).cartesian_begin(); x_it != complex_.point(v).cartesian_end(); ++x_it){
					if(first) first = false;
					else file<<",";
					file <<*x_it ;
				}
				file<<std::endl;
			}
		}
		else std::cerr << "Could not open "<<name_arff_file<<std::endl;

	}

	void read_points_from_arff_file(){
		std::ifstream file(name_arff_file_reduced);
		if(file.is_open()){
			std::string current_line;
			//goto @DATA
			do{
				std::getline(file,current_line);
			} while(current_line.find("@DATA") == std::string::npos);
			std::cout <<"points\n";
			while(std::getline(file,current_line)){
				if(current_line.find(",")!=std::string::npos){ //ignore lines without ,
					boost::char_separator<char> sep(",");
					boost::tokenizer<boost::char_separator<char>> tokens(current_line, sep);

					std::vector<double> point;
					for (const auto& t : tokens)
						point.push_back(std::stod(t));
					assert(point.size()==reduced_dimension_);
					// add a 0 coordinates to have at least 3 coordinates
					while(point.size()<3) point.push_back(0);
					reduced_points_.push_back(Point_3(point[0],point[1],point[2]));
					std::cout<<reduced_points_.back()<<std::endl;
				}
			}
		}
		else std::cerr << "Could not open "<<name_arff_file<<std::endl;
	}

	std::string name_arff_file = "tmp.arff";
	std::string name_arff_file_reduced = "tmp_red.arff";

	const Complex& complex_;
	std::vector<Point_3> reduced_points_;
	int reduced_dimension_;
};



#endif /* DIMENSION_REDUCTION_H_ */
