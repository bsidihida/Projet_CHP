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

#include "ILevelSet.hpp"
#include "Transport.hpp"
#include <limits>
#include "Utils.hpp"

/**
 * @file ILevelSet.cpp
 *
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */
namespace neos {

double ILevelSet::getLevelSet(const NPoint & cell) {
  double phi =  std::numeric_limits<double>::max();

  if (_geo.size() > 1) {
    for (size_t i = 0; i < _geo.size(); i++) {
      if (phi > _geo.at(i)->getLevelSet(cell)) {
        phi = _geo.at(i)->getLevelSet(cell);
      }
    }
  } else {
    _geo.front()->getLevelSet(cell);
  }
  return phi;
}

double ILevelSet::getLevelSet(const long &id) {
  double phi =  std::numeric_limits<double>::max();

  for (size_t i = 0; i < _geo.size(); i++) {
    double tphi = _geo[i]->getLevelSet(id);
    if (tphi < phi) {
      phi = tphi;
    }
  }
  return phi;
}

std::vector<double> ILevelSet::getLevelSet() {
  PiercedVector<double> fphi;

  fphi.reserve(_grid->nbCells());
  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior()) {
      const long &i=cell.getId();
      fphi.emplace(i);
    }
  }
  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior()) {
      const long &id = cell.getId();
      fphi[id] = getLevelSet(id);
    }
  }

  std::vector<double> phi = PVtoV(fphi, _grid);
  return phi;
}

void ILevelSet::transport(Grid *grid, const PiercedVector<NPoint> &u, const double dt) {
  Transport *trpt = new Transport(grid);

  for (size_t i = 0; i < _geo.size(); i++) {
    trpt->compute(_geo[i]->getPhi(), u, dt);
  }
  delete trpt;
}
}
