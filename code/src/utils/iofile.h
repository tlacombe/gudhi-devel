/*
 *  iofile.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/9/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */
#ifndef GUDHI_IOFILE_H_
#define GUDHI_IOFILE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <list>
#include <cassert>
#include <boost/graph/adjacency_list.hpp>


/**
 * \brief Read a set of points to turn it
 * into a vector< vector<double> > by filling points
 *
 * File format: 1 point per line
 * X11 X12 ... X1d 
 * X21 X22 ... X2d
 * etc
 */
inline void
read_points ( std::string file_name
            , std::vector< std::vector< double > > & points)
{  
  std::ifstream in_file (file_name.c_str(),std::ios::in);
  if(!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;}

  std::string line;
  double x;
  while( getline ( in_file , line ) )
  {
    std::vector< double > point;
    std::istringstream iss( line );
    while(iss >> x) { point.push_back(x); }
    points.push_back(point);
  }
  in_file.close();
}

/**
 * \brief Read a graph from a file.
 *
 * File format: 1 simplex per line
 * Dim1 X11 X12 ... X1d Fil1 
 * Dim2 X21 X22 ... X2d Fil2
 * etc
 *
 * The vertices must be labeled from 0 to n-1.
 * Every simplex must appear exactly once.
 * Simplices of dimension more than 1 are ignored.
 */
struct edge_filtration_tag {
  typedef boost::edge_property_tag kind;
};
struct vertex_filtration_tag {
  typedef boost::vertex_property_tag kind;
};

typedef int    Vertex_handle;
typedef double Filtration_value;
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
                              , boost::property < vertex_filtration_tag, Filtration_value >
                              , boost::property< edge_filtration_tag, Filtration_value > 
                              >                   Graph_t;
typedef std::pair< Vertex_handle, Vertex_handle > Edge_t;
inline Graph_t
read_graph ( std::string file_name )
{
  std::ifstream in_ (file_name.c_str(),std::ios::in);
  if(!in_.is_open()) { std::cerr << "Unable to open file " << file_name << std::endl; }

  std::vector< Edge_t >                       edges;
  std::vector< Filtration_value >             edges_fil;
  std::map< Vertex_handle, Filtration_value > vertices;
  
  std::string      line;
  int              dim;
  Vertex_handle    u,v,max_h = -1;
  Filtration_value fil;
  while( getline ( in_ , line ) )
  {
    std::istringstream iss( line );
    while(iss >> dim) {
      switch ( dim ) {
        case 0 : {
          iss >> u; iss >> fil; 
          vertices[u] = fil;
          if(max_h < u) { max_h = u; }
          break;
        }
        case 1 : {
          iss >> u; iss >> v; iss >> fil;
          edges.push_back(Edge_t(u,v));
          edges_fil.push_back(fil);
          break;
        }
        default: {break;} 
      } 
    }
  }
  in_.close();
  if(max_h+1 != vertices.size()) 
    { std::cerr << "Error: vertices must be labeled from 0 to n-1 \n"; }

  Graph_t skel_graph(edges.begin(),edges.end(),edges_fil.begin(),vertices.size());
  auto vertex_prop = boost::get(vertex_filtration_tag(),skel_graph);

  boost::graph_traits<Graph_t>::vertex_iterator vi, vi_end;
  auto v_it = vertices.begin();
  for (tie(vi, vi_end) = boost::vertices(skel_graph); vi != vi_end; ++vi,++v_it)
  { boost::put(vertex_prop,*vi,v_it->second); }

  return skel_graph;
}

/**
 * @brief OFF reader, save the content of the file (vertices and maximal faces) in 'complex'.
 * The class Complex has to handle the following operations:
 * - void add_vertex(double[dim],int dim)
 * - void add_face(int dimension ,int[dimension] vertices)
 * Source from : http://www.holmes3d.net/graphics/offfiles/OFFLoading.txt
 *
 * @todo todo adapt for points of arbitrary dimensions
 */
template<typename Complex> inline
bool general_read_off_file(const std::string & file_name, Complex& complex){
  // Declare temporary variables to read data into.
  // If the read goes well, we'll copy these into
  // our class variables, overwriting what used to
  // be there. If it doesn't, we won't have messed up
  // our previous data structures.
  int tempNumPoints   = 0;  // Number of x,y,z coordinate triples
  int tempNumFaces    = 0;  // Number of polygon sets
  int tempNumEdges    = 0;  // Unused, except for reading.
  int tempDimPoints   = 0;
  double** tempPoints = NULL;  // An array of x,y,z coordinates.
  int** tempFaces = NULL;    // An array of arrays of point
  // pointers. Each entry in this
  // is an array of integers. Each
  // integer in that array is the
  // index of the x, y, and z
  // coordinates in the corresponding
  // arrays.
  int*  tempFaceSizes = NULL;  // An array of polygon point counts.
  // Each of the arrays in the tempFaces
  // array may be of different lengths.
  // This array corresponds to that
  // array, and gives their lengths.
  int i;        // Generic loop variable.
  bool goodLoad = true;    // Set to false if the file appears
  // not to be a valid OFF file.
  char tempBuf[128];    // A buffer for reading strings
  // from the file.

  // Create an input file stream for the file the CArchive
  // is connected to. This allows use of the overloaded
  // extraction operator for conversion from string
  // data to doubles and ints.
  std::ifstream ifs (file_name.c_str(), std::ios::in);
  if(!ifs.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return false;
  }

  // Grab the first string. If it's "OFF", we think this
  // is an OFF file and continue. Otherwise we give up.
  ifs >> tempBuf;
  if (strcmp(tempBuf, "OFF") != 0) {
    goodLoad = false;
    std::cerr << "No OFF preambule\n";
  }

  // Read the sizes for our two arrays, and the third
  // int on the line. If the important two are zero
  // sized, this is a messed up OFF file. Otherwise,
  // we setup our temporary arrays.
  if (goodLoad) {
    ifs >> tempNumPoints >> tempNumFaces >> tempNumEdges;
    if (tempNumPoints < 1 || tempNumFaces < 0) {
      // If either of these were negative, we make
      // sure that both are set to zero. This is
      // important for later deleting our temporary
      // storage.
      goodLoad      = false;
      std::cerr << "tempNumPoints < 1 || tempNumFaces < 0\n";
      tempNumPoints = 0;
      tempNumFaces  = 0;
    } else {
      tempPoints = new double*[tempNumPoints];
      tempFaces = new int*[tempNumFaces];
      tempFaceSizes = new int[tempNumFaces];
    }
  }

  if (goodLoad) {
    // Load all of the points.

    // we start by loading the first one
    // the case is difference because then we dont know the dimension by advance
    // we discover the point dimension by reading the first line
    std::string lineFirstPoint;
    std::getline(ifs, lineFirstPoint);
    assert(lineFirstPoint.empty());
//    xxx do sthing else, incompatible with clang compil
//    goodLoad =
    std::getline(ifs, lineFirstPoint);
//    if(!goodLoad)
//      std::cerr<<"Cant read the first point\n";
    std::istringstream lineFirstPointStream(lineFirstPoint);

    // we store the first point in a temporary list
    std::list<double> firstTempPoint;
    double coord;
    while(lineFirstPointStream>>coord){
      firstTempPoint.push_back(coord);
      ++tempDimPoints;
    }
    // we store the point in our points array
    tempPoints[0]=new double[tempDimPoints];
    for( int j = 0 ; j<tempDimPoints; ++j){
      tempPoints[0][j] = firstTempPoint.front();
      firstTempPoint.pop_front();
    }

    // now, we know the dimension and can read safely other points
    for (i = 1; i < tempNumPoints; i++) {
      tempPoints[i] = new double[tempDimPoints];
      for (int j = 0 ; j<tempDimPoints ; ++j){
        ifs>>tempPoints[i][j];
      }
    }

    // Load all of the faces.
    for (i = 0; i < tempNumFaces; i++) {
      // This tells us how many points make up
      // this face.
      ifs >> tempFaceSizes[i];
      // So we declare a new array of that size
      tempFaces[i] = new int[tempFaceSizes[i]];
      // And load its elements with the vertex indices.
      for (int j = 0; j < tempFaceSizes[i]; j++) {
        ifs >> tempFaces[i][j];
      }
      // Clear out any face color data by reading up to
      // the newline. 128 is probably considerably more
      // space than necessary, but better safe than
      // sorry.
      ifs.getline(tempBuf, 128);
    }
  }

  // Here is where we copy the data from the temp
  // structures into our permanent structures. We
  // probably will do some more processing on the
  // data at the same time. This code you must fill
  // in on your own.
  if (goodLoad) {
    // we save vertices first in the complex
    for (i = 0; i < tempNumPoints; i++)
      complex.add_vertex(tempPoints[i],tempDimPoints);

    // we save faces
    for (i = 0; i < tempNumFaces; i++) {
      for (int j = 0; j < tempFaceSizes[i]; j++)
        complex.add_face(tempFaceSizes[i],tempFaces[i]);
    }
  }

  // Now that we're done, we have to make sure we
  // free our dynamic memory.
  for (i = 0; i < tempNumPoints; i++) {
    delete []tempPoints[i];
  }
  delete []tempPoints;

  for (i = 0; i < tempNumFaces; i++) {
    delete tempFaces[i];
  }
  delete []tempFaces;
  delete []tempFaceSizes;

  // Clean up our ifstream. The MFC framework will
  // take care of the CArchive.
  ifs.close();

  return goodLoad;
}


template<typename Complex>
class Geometric_flag_complex_wrapper{
  Complex& complex_;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Point Point;

  const bool load_only_points_;

public:
  Geometric_flag_complex_wrapper(Complex& complex,bool load_only_points = false):
    complex_(complex),
    load_only_points_(load_only_points)
{}


  void add_vertex(double* xyz,int dim){
    Point p(dim);
    for(int i=0;i<dim;++i)
      p[i] = xyz[i];
    complex_.add_vertex(p);
  }

  void add_face(int dimension ,int* vertices){
    if (!load_only_points_){
      for (int i = 0; i<dimension ; ++i)
        for (int j = i+1; j<dimension ; ++j)
          complex_.add_edge(Vertex_handle(vertices[i]),Vertex_handle(vertices[j]));
    }
  }
};




/**
 * @brief Read a mesh into a OFF file
 * load_only_points should be true if only the points have to be loaded.
 */
template<typename Complex>
bool read_off_file(std::string file_name,Complex &complex,bool load_only_points = false){
  complex.clear();
  Geometric_flag_complex_wrapper<Complex> complex_wrapper(complex,load_only_points);
  return general_read_off_file(file_name,complex_wrapper);

}


#endif
