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
 * @file Geo2D.cpp
 * @brief This file contains 2D Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-20
 * @copyright Inria
 */

#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Geo2D.hpp"

namespace neos {

void Geo2D::levelSet(Grid *grid) {
  std::vector<double> tau;

  _phi.clear();
  _phi.reserve(grid->getCellCount());
  for (auto cell = grid->internalCellBegin(); cell != grid->internalCellEnd();
       ++cell) {
    const long &i=cell.getId();
    _phi.emplace(i);
  }

  for (int k = 0; k < _precision; k++) {
    double taux = 0;
    double tauy = 0;

    taux = _bdr[k + 1][0] - _bdr[k][0];
    tauy = _bdr[k + 1][1] - _bdr[k][1];
    tau.push_back(taux/sqrt(pow(taux, 2) + pow(tauy, 2)));
    tau.push_back(tauy/sqrt(pow(taux, 2) + pow(tauy, 2)));
  }

  for (auto cell = grid->internalCellBegin(); cell != grid->internalCellEnd();
       ++cell) {
    double dist = 100000;
    const long &id = cell.getId();
    std::array<double,3> x0 = grid->evalCellCentroid(id);
    double px = 0;
    double py = 0;
    int kmin = 0;

    for (int k = 0; k < _precision; k++) {
      int xc = (_bdr[k + 1][0] - _bdr[k][0]) / 2;
      int yc = (_bdr[k + 1][1] - _bdr[k][1]) / 2;
      int distnew = sqrt(pow(x0[0] - xc, 2) + pow(x0[1] - yc, 2));

      if (distnew <= dist) {
        dist = distnew;
        px = x0[0] - xc;
        py = x0[1] - yc;
        kmin = k;
      }
    }
    _phi[id] = (tau[kmin] * py - tau[kmin+1] * px) / fabs(tau[kmin] * py - tau[kmin+1] * px) * dist;
  }
}

double Geo2D::getLevelSet(const std::array<double, 3> &cell) {
  std::vector<double> tau;

  for (int k = 0; k < _precision; k++) {
    double taux = 0;
    double tauy = 0;

    taux = _bdr[k + 1][0] - _bdr[k][0];
    tauy = _bdr[k + 1][1] - _bdr[k][1];
    tau.push_back(taux/sqrt(pow(taux, 2) + pow(tauy, 2)));
    tau.push_back(tauy/sqrt(pow(taux, 2) + pow(tauy, 2)));
  }

  double dist = 100000;
  double px = 0;
  double py = 0;
  int kmin = 0;

  for (int k = 0; k < _precision; k++) {
    int xc = (_bdr[k + 1][0] - _bdr[k][0]) / 2;
    int yc = (_bdr[k + 1][1] - _bdr[k][1]) / 2;
    int distnew = sqrt(pow(cell[0] - xc, 2) + pow(cell[1] - yc, 2));

    if (distnew <= dist) {
      dist = distnew;
      px = cell[0] - xc;
      py = cell[1] - yc;
      kmin = k;
    }
  }
  return (tau[kmin] * py - tau[kmin+1] * px) / fabs(tau[kmin] * py - tau[kmin+1] * px) * dist;
}
}
