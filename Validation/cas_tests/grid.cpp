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

#include "Grid.hpp"
#include "gtest/gtest.h"

using namespace neos;

TEST(Grid, init)
{
  Grid *grid = new Grid(0, 0, 0, 1.0, 0.5, 2);

  int nbCells = grid->nbCells();
  EXPECT_EQ(4, nbCells);

  int dim = grid->getDim();
  EXPECT_EQ(2, dim);

  delete grid;
}


TEST(Grid, locate)
{
  Grid *grid = new Grid(0, 0, 0, 3, 1., 2);
  grid->write();
  long idPoint = grid->locatePoint({ {  0, 0, 0} } );
  EXPECT_EQ(0, idPoint);

  std::vector<long> list = grid->findQuadrantNeighs({ { 1, 1, 0} }, 0);
  EXPECT_EQ(0, list[0]);
  EXPECT_EQ(1, list[1]);
  EXPECT_EQ(2, list[2]);
  EXPECT_EQ(3, list[3]);
  delete grid;
}
