/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Primoz Skraba
 *
 *    Copyright (C) 2009 Primoz Skraba.  All Rights Reserved.
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

#ifndef SRC_TOMATO_INCLUDE_GUDHI_TOMATO_INTERVAL_H_
#define SRC_TOMATO_INCLUDE_GUDHI_TOMATO_INTERVAL_H_

#include <cassert>

//======================
// basic data structure
// for holding birth and
// death times
//======================

class Interval {
 private:
  double birth;
  double death;
  bool infinite;

 public:
  Interval() { }

  Interval(double birth_) {
    birth = birth_;
    infinite = true;
  }

  void close(double death_) {
    death = death_;
    infinite = false;
  }

  double getBirth() {
    return birth;
  }

  double getDeath() {
    assert(!infinite);
    return death;
  }

  bool isInfinite() {
    return infinite;
  }
};

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_INTERVAL_H_
