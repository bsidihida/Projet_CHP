/* -------------------------------------------------------------------------*\
 *
 *  NEOS
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of Neos.
 *
 *  Neos is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  Neos is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Neos. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

/**
 * @file Sphere.cpp
 * @brief This file contains the Circle class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <cmath>
#include "Sphere.hpp"

namespace neos {

void Sphere::compute() {
  double pi = 2 * acos(0.);
  double phi = 0;
  double dphi = (2 * pi) / _precision;
  double theta = 0;
  double dtheta = pi / _slice;
  double radius = _radius;

  for (int k = 0; k < _precision; k++) {
    for (int i = 0; i < _slice; i++) {
      std::array<double, 3> pos;

      pos[0] = radius * sin(theta)  * cos(phi) + _center[0];
      pos[1] = radius * sin(theta)  * sin(phi) + _center[1];
      pos[2] = radius * cos(theta) + _center[2];
      _bdr.push_back(pos);
      theta += dtheta;
    }
    phi += dphi;
  }
}
}
