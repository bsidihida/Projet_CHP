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

#include <fstream>
#include "STLGeometry.hpp"
#include "surfunstructured.hpp"
#include "levelSetObject.hpp"

namespace neos {

STLGeometry::STLGeometry(const std::string &file, Grid *grid) {
  addSTL(file, grid);
  updatePhi();
}

STLGeometry::~STLGeometry() {
}

double STLGeometry::getLevelSet(const std::array<double, 3> &point) {
  return _phi[_grid->locatePoint(point)];
}

void STLGeometry::computeLevelSet(Grid *grid) {
  _grid = grid;
  _grid->update(true);
  _lvst.update(_grid->getMapper());
  updatePhi();
}

void STLGeometry::addSTL(const std::string &file, Grid *grid) {
  _grid = grid;
  bitpit::SurfUnstructured _stl(-1,3,3,MPI_COMM_WORLD);

  std::ifstream ifile(file);

  if (!ifile) {
    std::cout << "Cannot find : " << file << std::endl;
    exit(1);
  }

  if (_stl.importSTL(file, true) > 0 && _stl.importSTL(file, false) > 0) {
    std::cout << "import error" << std::endl;
    exit(1);
  }

  _lvst.setMesh(_grid);
  _stlID = _lvst.addObject(&_stl, 0);
  _lvst.compute();
}

void STLGeometry::updatePhi() {
  _phi.clear();
  _phi.reserve(_grid->nbCells());
  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior()) {
      const long &i=cell.getId();
      _phi.emplace(i);
    }
  }
  const bitpit::LevelSetObject &object = _lvst.getObject(_stlID);
  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior()) {
      const long &id = cell.getId();
      _phi[id] = object.getLS(id);
    }
  }
}
}
