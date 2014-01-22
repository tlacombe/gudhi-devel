/*
 *  iofile.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/9/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * \brief Read a set of points to turn it
 * into a vector< vector<double> > by filling points
 *
 * File format: 1 point per line
 * X11 X12 ... X1d 
 * X21 X22 ... X2d
 * etc
 */
void
read_points(std::string file_name,
	    std::vector< std::vector< double > > &points)
{	
  std::ifstream in_file (file_name.c_str(),std::ios::in);
  if(!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;}
	
  std::string line;
  double x;
  while( getline (in_file,line) )
    {
      std::vector<double> point;
      std::istringstream iss(line);
      while(iss >> x) {point.push_back(x);}
      points.push_back(point);		
    }
  in_file.close();
}
