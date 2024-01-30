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
 * @file ACylinder.cpp
 * @brief This file contains Analytic Cylinder class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <iostream>
#include <fstream>
#include "ACylinder.hpp"

namespace neos {

double ACylinder::getLevelSet(const NPoint &cell) {
  std::vector<double> v;

  v.resize(_dim, 0.0);
  for (int j = 0; j < _dim; ++j) {
    if (j != _axe) {
      v[j] = cell[j] - _center[j];
    }
  }
  return norm2(v) - _radius;
}

void ACylinder::computeLevelSet(Grid *g) {
  NPoint cellcenter;
  std::vector<double> v;

  v.resize(_dim, 0.0);
  _phi.clear();
  _phi.reserve(g->getCellCount());
  for (auto cell = g->internalCellBegin(); cell != g->internalCellEnd(); ++cell) {
    const long &i=cell.getId();
    _phi.emplace(i);
  }

  for (auto cell = g->internalCellBegin(); cell != g->internalCellEnd(); ++cell) {
    const long &id = cell.getId();
    cellcenter = g->evalCellCentroid(id);
    for (int j = 0; j < _dim; ++j) {
      if (j != _axe) {
        v[j] = cellcenter[j] - _center[j];
      }
    }
    _phi[id] = norm2(v) - _radius;
  }
}
}
