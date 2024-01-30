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
 * @file Geo3D.cpp
 * @brief This file contains 3D Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-20
 * @copyright Inria
 */

#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Geo3D.hpp"
#include "Grid.hpp"

namespace neos {

void Geo3D::levelSet(Grid *grid) {
  _phi.clear();
  _phi.reserve(grid->getCellCount());
  for (auto cell = grid->internalCellBegin(); cell != grid->internalCellEnd();
       ++cell) {
    const long &i=cell.getId();
    _phi.emplace(i);
  }

  for (auto cell = grid->internalCellBegin(); cell != grid->internalCellEnd();
       ++cell) {
    const long &id = cell.getId();
    std::array<double,3> x0 = grid->evalCellCentroid(id);
    double dist = 100000;

    for (int n = 0; n < _slice - 1; n++) {
      for (int l = 0; l < _precision - 1; l++) {
        int xc = (_bdr[(_slice * n) + l][0] + _bdr[(_slice * (n+1)) + l][0]  +
                  _bdr[(_slice * n) + l + 1][0] + _bdr[(_slice * (n+1)) + l + 1][0]) / 4;
        int yc = (_bdr[(_slice * n) + l][1] + _bdr[(_slice * (n+1)) + l][1]  +
                  _bdr[(_slice * n) + l + 1][1] + _bdr[(_slice * (n+1)) + l + 1][1]) / 4;
        int zc = (_bdr[(_slice * n) + l][2] + _bdr[(_slice * (n+1)) + l][2] +
                  _bdr[(_slice * n) + l + 1][2] + _bdr[(_slice * (n+1)) + l + 1][2]) / 4;

        int distnew = sqrt(pow(x0[0] - xc, 2) + pow(x0[1] - yc, 2) + pow(x0[2] - zc, 2));

        if (distnew <= dist)
          dist = distnew;
      }
    }
    _phi[id] = dist;
  }
}

double Geo3D::getLevelSet(const std::array<double, 3> &cell) {
  double dist = 100000;

  for (int n = 0; n < _slice - 1; n++) {
    for (int l = 0; l < _precision - 1; l++) {
      int xc = (_bdr[(_slice * n) + l][0] + _bdr[(_slice * (n+1)) + l][0]  +
                _bdr[(_slice * n) + l + 1][0] + _bdr[(_slice * (n+1)) + l + 1][0]) / 4;
      int yc = (_bdr[(_slice * n) + l][1] + _bdr[(_slice * (n+1)) + l][1]  +
                _bdr[(_slice * n) + l + 1][1] + _bdr[(_slice * (n+1)) + l + 1][1]) / 4;
      int zc = (_bdr[(_slice * n) + l][2] + _bdr[(_slice * (n+1)) + l][2] +
                _bdr[(_slice * n) + l + 1][2] + _bdr[(_slice * (n+1)) + l + 1][2]) / 4;

      int distnew = sqrt(pow(cell[0] - xc, 2) + pow(cell[1] - yc, 2) + pow(cell[2] - zc, 2));

      if (distnew <= dist)
        dist = distnew;
    }
  }
  return dist;
}
}
