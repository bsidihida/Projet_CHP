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

#include "Neos.hpp"
#include "gtest/gtest.h"
#include <math.h>

using namespace neos;

TEST(LevelSet, addGeo) {
  Grid        *grid = new Grid(0, 0, 0, 1.0, 0.1, 3);
  ASphere     *geo  = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo2 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo3 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  int idGeo2;

  LevelSet ls(grid);
  ls.addGeometry(geo);
  idGeo2 = ls.addGeometry(geo2);
  ls.addGeometry(geo3);

  Geometry *res = ls.getGeometry(idGeo2);

  EXPECT_EQ(res, geo2);

  delete geo;
  delete geo2;
  delete geo3;
  delete grid;
}

TEST(LevelSet, addGeo2) {
  Grid        *grid = new Grid(0, 0, 0, 1.0, 0.1, 3);
  ASphere     *geo  = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo2 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo3 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  int idGeo2;
  int idGeo3;

  LevelSet ls(grid);
  ls.addGeometry(geo);
  idGeo2 = ls.addGeometry(geo2, "test");
  idGeo3 = ls.addGeometry(geo3, "test");

  EXPECT_EQ(idGeo2, idGeo3);
  EXPECT_EQ(ls.countGeometry(), 2);

  int res = ls.getId("test");

  EXPECT_EQ(idGeo2, res);

  delete geo;
  delete geo2;
  delete geo3;
  delete grid;
}

TEST(LevelSet, delGeo) {
  Grid        *grid = new Grid(0, 0, 0, 1.0, 0.1, 3);
  ASphere     *geo  = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo2 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo3 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  ASphere     *geo4 = new ASphere(0.5, 0.5, 0.0, 0.2, 3);
  int idGeo;
  int idGeo2;

  LevelSet ls(grid);
  idGeo  = ls.addGeometry(geo);
  ls.addGeometry(geo2, "test");
  ls.addGeometry(geo4, "test");
  ls.addGeometry(geo3, "test");
  EXPECT_EQ(2, ls.countGeometry());

  ls.delGeometry("test");
  EXPECT_EQ(1, ls.countGeometry());

  ls.delGeometry(idGeo);
  EXPECT_EQ(0, ls.countGeometry());

  idGeo  = ls.addGeometry(geo);
  idGeo2 = ls.addGeometry(geo2, "test");

  ls.delGeometry(idGeo);
  EXPECT_EQ(1, ls.countGeometry());

  ls.delGeometry(idGeo2);
  EXPECT_EQ(0, ls.countGeometry());
  delete geo;
  delete geo2;
  delete geo3;
  delete grid;
}

TEST(LevelSet, ASphere) {
  Grid        *grid = new Grid(0, 0, 0, 1.0, 0.1, 3);
  ASphere     *geo  = new ASphere(0.5, 0.5, 0.5, 0.2, 3);

  LevelSet ls(grid);
  ls.addGeometry(geo);

  geo->computeLevelSet(grid);

  std::vector<double> v(grid->getCellCount());

  std::vector<double> phiV = ls.getLevelSet();

  bitpit::PiercedVector<double> phi;
  NPoint cellcenter;

  phi.reserve(grid->getCellCount());
  for (auto &cell : grid->getCells()) {
    const long &i = cell.getId();
    phi.emplace(i);
  }

  for (auto &cell : grid->getCells()) {
    const long &id = cell.getId();
    cellcenter = grid->evalCellCentroid(id);
    for (int j = 0; j < 3; ++j) {
      v[j] = cellcenter[j] - 0.5;
    }
    phi[id] = norm2(v) - 0.2;
  }

  double error = 0;
  double error2 = 0;
  double error3 = 0;
  for (auto &cell : grid->getCells()) {
    const long &id = cell.getId();
    error = std::max((double)abs(phi[id] - phiV[id]), error);
    error2 = std::max((double)abs(phi[id] - geo->getLevelSet(grid->evalCellCentroid(id))), error2);
    error3 = std::max((double)abs(phi[id] - geo->getLevelSet(id)), error3);
  }
  EXPECT_EQ(0, error);
  EXPECT_EQ(0, error2);
  EXPECT_EQ(0, error3);
}

TEST(LevelSet, ACylinder) {
  Grid        *grid = new Grid(0, 0, 0, 1.0, 0.1, 3);
  ACylinder   *geo  = new ACylinder(0.5, 0.5, 0.5, 0.2, 3, 2);

  LevelSet ls(grid);
  ls.addGeometry(geo);

  geo->computeLevelSet(grid);

  std::vector<double> v(grid->getCellCount());

  std::vector<double> phiV = ls.getLevelSet();

  bitpit::PiercedVector<double> phi;
  NPoint cellcenter;

  phi.reserve(grid->getCellCount());
  for (auto &cell : grid->getCells()) {
    const long &i = cell.getId();
    phi.emplace(i);
  }

  for (auto &cell : grid->getCells()) {
    const long &id = cell.getId();
    cellcenter = grid->evalCellCentroid(id);
    for (int j = 0; j < 3; ++j) {
      if (j != 2) {
        v[j] = cellcenter[j] - 0.5;
      }
    }
    phi[id] = norm2(v) - 0.2;
  }

  double error = 0;
  double error2 = 0;
  double error3 = 0;
  for (auto &cell : grid->getCells()) {
    const long &id = cell.getId();
    error = std::max((double)abs(phi[id] - phiV[id]), error);
    error2 = std::max((double)abs(phi[id] - geo->getLevelSet(grid->evalCellCentroid(id))), error2);
    error3 = std::max((double)abs(phi[id] - geo->getLevelSet(id)), error3);
  }
  EXPECT_EQ(0, error);
  EXPECT_EQ(0, error2);
  EXPECT_EQ(0, error3);
}
