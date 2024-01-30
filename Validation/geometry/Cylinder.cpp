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
 * @file Cylinder.cpp
 * @brief This file contains the Cylinder class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <cmath>
#include "Cylinder.hpp"

namespace neos {

void Cylinder::compute() {
  double beta = 0;
  double dbeta = 0;
  double pi = 2 * acos(0.);
  double r  = _radius;

  dbeta = (2*pi) / _precision;
  for(int k = 0; k < _slice; k++) {
    for (int l = 0; l < _precision; l++) {
      std::array<double, 3> pos;

      pos[0] = r*cos(beta) + _center[0];
      pos[1] = r*sin(beta) + _center[1];
      pos[2] = k*(_lenght/_slice) + _center[2];
      _bdr.push_back(pos);
      beta += dbeta;
    }
  }
}
}
