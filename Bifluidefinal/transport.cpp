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
#include <fstream>
#include <vector>
#include <iterator>
#include <array>
#include <math.h>
#include <string>
#include <algorithm>
#include "Neos.hpp"

/**
 * Transport.cpp
 * Example to show how to use neos transport class
 **/

/*
 * list of command line arguments
 * H <- cell size
 * ITYPE <- Interpolation type (see interpo_type enum)
 * STEPS <- Number of steps
 * DIM <- Grid dimension
 * EPSILON <- For RBF Interpolation (optional)
 *
 * NBAND <- Narrow Band size for error calculation step
 */
#define H 1
#define ITYPE 2
#define STEPS 3
#define DIM 4
#define EPSILON 5
#define NBAND 0.1

using namespace neos;

int main(int ac, char **av) {
  Neos_Init();

  double h        = std::stod(av[H]);
  int itype    = std::stoi(av[ITYPE]);
  int nbSteps  = std::stoi(av[STEPS]);
  int dim      = std::stoi(av[DIM]);
  // Neos analytic sphere
  ASphere   *geo      = new ASphere(0.5, 0.5, 0, 0.2, dim);
  // Neos grid
  Grid      *grid     = new Grid(0, 0, 0, 1.0, h, dim);
  double simuStep = h / 0.3 * 0.8;
  // Neos transport class
  Transport *trpt     = new Transport(grid);

  // If RBF it's used, you can set an epsilon value
  if (itype == RBF && ac == EPSILON) {
    trpt->setEpsilon(std::stod(av[EPSILON]));
  }

  // Speed vector
  std::vector<double> u;
  u.push_back(0.3);
  u.push_back(0.3);
  u.push_back(0.3);

  int nocts = grid->nbCells();
  grid->setExportName("transport");

  /**
   * Compute the levelset
   * make a step
   * update the geometry center to calculate the error
   **/
  geo->computeLevelSet(grid);
  for (int i = 0; i < nbSteps; i++) {
    trpt->compute(geo->getPhi(), u, simuStep);
    geo->updateCenter(simuStep, u);
  }

  std::vector<double> ephi(nocts,0.0);
  for (auto & cell : grid->getCells()) {
    const long &id = cell.getId();
    ephi[id] = geo->getLevelSet(id);
  }

  double error = 0;
  std::vector<double> eephi(nocts,0.0);
  std::vector<double>::iterator eitphi = eephi.begin();
  std::vector<double> phi(nocts, 0.0);
  std::vector<double> aphi(nocts,0.0);
  std::vector<double>::iterator aitphi = aphi.begin();
  std::vector<double>::iterator itphi = phi.begin();
// compute the levelset
  geo->computeLevelSet(grid);

  // copy the levelset value to a vector
  for (auto & cell : grid->getCells()) {
    const long &id = cell.getId();
    *eitphi = abs(ephi[id] - geo->getLevelSet(id));
    *itphi = geo->getLevelSet(id);

    // Use the narrow band for error calculation
    if (abs(ephi[id]) <= NBAND) {
      error = std::max((double)abs(ephi[id] - geo->getLevelSet(id)), error);
    }
    ++eitphi;
    *aitphi = geo->getLevelSet(id);
    ++aitphi;
    ++itphi;
  }
// Export to file
  phi = PVtoV(geo->getPhi(),grid);
  grid->addData("levelset", phi);
  grid->addData("error", eephi);
  grid->addData("aphi", aphi);
  grid->write();

  delete geo;
  delete grid;
  delete trpt;

  Neos_Finalize();

  return 0;
}
