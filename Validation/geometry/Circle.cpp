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
 * @file Circle.cpp
 * @brief This file contains the Circle class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */


#include <cmath>
#include "Circle.hpp"

namespace neos {

void Circle::compute() {
  double pi = 2 * acos(0.);
  double beta = 0;
  double dbeta = 2 * pi / (_precision - 1);
  double radius = _radius;

  for (int i = 0; i < _precision; i++) {
    NPoint pos;

    pos[NPX] = radius * cos(beta) + _center[NPX];
    pos[NPY] = radius * sin(beta) + _center[NPY];
    pos[NPZ] = 0.0;
    _bdr.push_back(pos);
    beta += dbeta;
  }
}
}
