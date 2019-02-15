/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "read_perseus_file"
#include <boost/test/unit_test.hpp>

#include <gudhi/Bitmap_cubical_complex/read_perseus_file.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <cstdio>

BOOST_AUTO_TEST_CASE(read_perseus_unexisting_file) {
  std::vector<unsigned> sizes;
  std::vector<double> filtrations;
  std::vector<bool> directions;


  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("some_impossible_weird_file_name",
                                                                             sizes, filtrations, directions),
                     std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(read_perseus_inconsistent_file) {
  // Data to be read
  std::vector<unsigned> sizes;
  std::vector<double> filtrations;
  std::vector<bool> directions;

  remove("perseus.txt");
  std::ofstream perseus_file;
  perseus_file.open ("perseus.txt");
  perseus_file << "# Just comments and empty line\n";
  perseus_file << "\n";
  perseus_file << "# Another comment\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::invalid_argument);

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Data dimension without dimensions\n";
  perseus_file << "5\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::invalid_argument);

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Invalid data dimension\n";
  perseus_file << "abc\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::ios_base::failure);

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Data dimension\n";
  perseus_file << "3\n";
  perseus_file << "# Requires 3 dimensions, provide only 2\n";
  perseus_file << "10\n";
  perseus_file << "10\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::invalid_argument);

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Data dimension\n";
  perseus_file << "3\n";
  perseus_file << "# Requires 3 dimensions, provide an invalid one\n";
  perseus_file << "10\n";
  perseus_file << "# Invalid dimension\n";
  perseus_file << "abc\n";
  perseus_file << "10\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::ios_base::failure);

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Data dimension\n";
  perseus_file << "2\n";
  perseus_file << "# Requires 2 dimensions - provided\n";
  perseus_file << "3\n";
  perseus_file << "3\n";
  perseus_file << "# Requires 9 filtrations - only 5 provided\n";
  perseus_file << "1.\n";
  perseus_file << "2.\n";
  perseus_file << "3.\n";
  perseus_file << "4.\n";
  perseus_file << "5.\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::invalid_argument);

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Data dimension\n";
  perseus_file << "2\n";
  perseus_file << "# Requires 2 dimensions - provided\n";
  perseus_file << "3\n";
  perseus_file << "3\n";
  perseus_file << "# Requires 9 filtrations - provide an invalid one\n";
  perseus_file << "1.\n";
  perseus_file << "2.\n";
  perseus_file << "3.\n";
  perseus_file << "4.\n";
  perseus_file << "5.\n";
  perseus_file << "6.\n";
  perseus_file << "# Invalid filtration\n";
  perseus_file << "G\n";
  perseus_file << "8.\n";
  perseus_file << "9.\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::ios_base::failure);

}

BOOST_AUTO_TEST_CASE(perseus_file_read_inf_value) {
  // Data to be read
  std::vector<unsigned> sizes;
  std::vector<double> filtrations;
  std::vector<bool> directions;

  Gudhi::cubical_complex::read_perseus_style_file<double>("sinusoid.txt", sizes, filtrations, directions);

  std::cout << "sinusoid.txt has " << filtrations.size() << " values." << std::endl;
  BOOST_CHECK(filtrations.size() == 75);
  std::cout << "First value of sinusoid.txt is " << filtrations[0] << std::endl;
  BOOST_CHECK(filtrations[0] == 10.);
  // Next value
  std::cout << "Second value of sinusoid.txt is " << filtrations[1] << std::endl;
  BOOST_CHECK(filtrations[1] == std::numeric_limits<double>::infinity());

  std::cout << "sinusoid.txt has " << sizes.size() << " dimension." << std::endl;
  BOOST_CHECK(sizes.size() == 1);
  std::cout << "First dimension value of sinusoid.txt is " << sizes[0] << std::endl;
  BOOST_CHECK(sizes[0] == 75);

  std::cout << "sinusoid.txt has " << directions.size() << " direction." << std::endl;
  BOOST_CHECK(sizes.size() == 1);
  std::cout << "First direction value of sinusoid.txt is " << directions[0] << std::endl;
  BOOST_CHECK(directions[0] == false);

}

BOOST_AUTO_TEST_CASE(perseus_file_read_periodic_value) {
  // Data to be read
  std::vector<unsigned> sizes;
  std::vector<double> filtrations;
  std::vector<bool> directions;

  Gudhi::cubical_complex::read_perseus_style_file<double>("periodiccubicalcomplexdoc.txt", sizes, filtrations,
                                                          directions);

  std::cout << "periodiccubicalcomplexdoc.txt has " << sizes.size() << " dimension." << std::endl;
  BOOST_CHECK(sizes.size() == 2);
  std::cout << "First dimension value of periodiccubicalcomplexdoc.txt is " << sizes[0] << std::endl;
  BOOST_CHECK(sizes[0] == 3);
  std::cout << "Second dimension value of periodiccubicalcomplexdoc.txt is " << sizes[1] << std::endl;
  BOOST_CHECK(sizes[1] == 3);

  std::cout << "periodiccubicalcomplexdoc.txt has " << directions.size() << " direction." << std::endl;
  BOOST_CHECK(directions.size() == 2);
  std::cout << "First direction value of periodiccubicalcomplexdoc.txt is " << directions[0] << std::endl;
  BOOST_CHECK(directions[0] == true);
  std::cout << "Second direction value of periodiccubicalcomplexdoc.txt is " << directions[1] << std::endl;
  BOOST_CHECK(directions[1] == false);

  std::cout << "periodiccubicalcomplexdoc.txt has " << filtrations.size() << " values." << std::endl;
  BOOST_CHECK(filtrations.size() == (sizes[0] * sizes[1]));
}

BOOST_AUTO_TEST_CASE(read_perseus_too_many_filtrations_file) {
  // Data to be read
  std::vector<unsigned> sizes;
  std::vector<double> filtrations;
  std::vector<bool> directions;

  remove("perseus.txt");
  std::ofstream perseus_file;

  remove("perseus.txt");
  perseus_file.open ("perseus.txt");
  perseus_file << "# Data dimension\n";
  perseus_file << "2\n";
  perseus_file << "# Requires 2 dimensions - provided\n";
  perseus_file << "3\n";
  perseus_file << "3\n";
  perseus_file << "# Requires 9 filtrations - provide an invalid one\n";
  perseus_file << "1.\n";
  perseus_file << "2.\n";
  perseus_file << "3.\n";
  perseus_file << "4.\n";
  perseus_file << "5.\n";
  perseus_file << "6.\n";
  perseus_file << "7.\n";
  perseus_file << "8.\n";
  perseus_file << "9.\n";
  perseus_file << "# Too many data\n";
  perseus_file << "10.\n";
  perseus_file.close();

  BOOST_CHECK_THROW (Gudhi::cubical_complex::read_perseus_style_file<double>("perseus.txt",
                                                                             sizes, filtrations, directions),
                     std::ios_base::failure);

}
