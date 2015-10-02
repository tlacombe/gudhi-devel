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

#ifndef _CLUSTER_INTERVAL_H
#define _CLUSTER_INTERVAL_H

#include <cassert>

//======================
// basic data structure
// for holding birth and
// death times
//======================

namespace Gudhi {

namespace cluster {

/**
 * \brief Cluster_interval is a basic data structure for holding birth and death times.
 *
 */
class Cluster_interval {
 private:
  /** \brief Birth of the cluster interval.*/
  double birth_;
  /** \brief Death of the cluster interval.*/
  double death_;
  /** \brief Infinite status of the cluster interval (true if not dead, false otherwise).*/
  bool infinite_;

  // Private default constructor for not to be used.
  Cluster_interval() { }

 public:
  /** \brief Cluster_interval construction from its birth value.*/
  Cluster_interval(double birth) {
    birth_ = birth;
    infinite_ = true;
  }

  /** \brief Cluster_interval closure from its death value.*/
  void close(double death) {
    death_ = death;
    infinite_ = false;
  }

  /** \brief Cluster_interval birth value getter.*/
  double get_birth() const {
    return birth_;
  }

  /** \brief Cluster_interval death value getter.
   * \warning getter asserts if the cluster interval is infinite.
   */
  double get_death() const {
    assert(!infinite_);
    return death_;
  }

  /** \brief Cluster_interval infinite status getter.*/
  bool is_infinite() const {
    return infinite_;
  }
};

}  // namespace cluster

}  // namespace Gudhi

#endif  // _CLUSTER_INTERVAL_H
