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

#include "bitpit.hpp"
#include "Transport.hpp"
#include "Grid.hpp"
#include "gtest/gtest.h"
#include <math.h>
#include "Neos.hpp"

using namespace neos;

TEST(Transport, type) {
  Grid *grid = new Grid(0, 0, 0, 1.0, 0.2, 3);
  Transport *trpt = new Transport(grid);

  trpt->setInterpType(POLY3D);
  int interp = trpt->getInterpType();

  delete grid;
  delete trpt;
  ASSERT_EQ(interp, POLY3D);
}

TEST(Transport, epsilon) {
  Grid *grid = new Grid(0, 0, 0, 1.0, 0.2, 3);
  Transport *trpt = new Transport(grid);

  trpt->setEpsilon(3);
  double epsilon = trpt->getEpsilon();

  delete grid;
  delete trpt;
  ASSERT_EQ(epsilon, 3);
}

TEST(Transport, compute) {
  Grid *grid = new Grid(0, 0, 0, 10.0, 1., 3);
  Transport *trpt = new Transport(grid);
  ASphere *geo = new ASphere(5., 5., 5., 2., 3);
  std::vector<double> phi;

  grid->setExportName("testexport");

  std::vector<double> u;
  u.push_back(.3);
  u.push_back(.3);
  u.push_back(.3);

  geo->computeLevelSet(grid);
  for (int i = 0; i < 3; i++) {
    phi = PVtoV(geo->getPhi(),grid);
    grid->addData("levelset", phi);
    trpt->compute(geo->getPhi(), u, 0.8);
    geo->updateCenter(0.8, u);
  }

  phi = PVtoV(geo->getPhi(),grid);
  grid->addData("levelset", phi);
  grid->write();
  delete trpt;
  delete grid;
  delete geo;
}
