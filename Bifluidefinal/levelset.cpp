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

#include <iostream>
#include <string>
#include "Neos.hpp"

using namespace neos;

int main() {

  Neos_Init();

  int dim  = GRID_3D;
  Grid        *grid = new Grid(-20, -10, -5, 20.0, 0.5, dim);
  ASphere     *geo  = new ASphere(10, 15, 10, 2.0, dim);
  ASphere     *geo2 = new ASphere(5, 7, 13, 2.0, dim);
  ASphere     *geo3 = new ASphere(7, 7, 3, 2.0, dim);
  STLGeometry *geo4 = new STLGeometry();
  geo4->addSTL("./Dragonite.stl", grid);
  LevelSet *ls = new LevelSet(grid);

  geo->computeLevelSet(grid);
  geo2->computeLevelSet(grid);
  geo3->computeLevelSet(grid);
  geo4->computeLevelSet(grid);

  ls->addGeometry(geo, "TAG");
  ls->addGeometry(geo2, "TAG");
  ls->addGeometry(geo4, "TAG");
  ls->addGeometry(geo3);


  std::vector<double> phi = ls->getLevelSet("TAG");

  grid->setExportName("levelset");
  grid->addData("phi", phi);
  grid->write();
  delete geo;
  delete geo2;
  delete geo3;
  delete grid;

  Neos_Finalize();

  return 0;
}
