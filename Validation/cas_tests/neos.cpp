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

#include "NeosInterface.hpp"
#include "gtest/gtest.h"

using namespace neos;

TEST(NeosInterface, api) {
  NeosInterface *neos = new NeosInterface();
  neos->createGrid({ {0., 0., 0.} }, 10, 0.2);
  neos->addSphere({ {5., 5., 5.} }, 2.0, "geo");
  neos->write();
  delete neos;
}

TEST(NeosInterface, api2) {
  NeosInterface *neos = new NeosInterface();
  neos->createGrid({ {0., 0., 0.} }, 10, 0.2);
  neos->addSphere({ {5., 5., 5.} }, 2.0, "geo");
  neos->write();

  std::vector<double> u =  { { 0.3, -0.3, 0.3} };

  for (int i = 0; i < 2; i++) {
    neos->transport("geo", u, .5 / .3 * .8);
    neos->write();
  }

  delete neos;
}


TEST(NeosInterface, api3) {
  NeosInterface *neos = new NeosInterface();
  neos->createGrid({ {0., 0., 0.} }, 10, 0.2);
  neos->addSphere({ {5., 5., 5.} }, 2.0, "geo");
  neos->write();

  std::vector<NPoint> u = neos->randomNPointVector();

  for (int i = 0; i < 2; i++) {
    neos->transport("geo", u, .5 / .3 * .8);
    neos->write();
  }

  delete neos;
}
