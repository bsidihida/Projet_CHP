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
#include "LaplacianFactory.hpp"
#include "Laplacian.hpp"
#include "Grid.hpp"
#include "Utils.hpp"
#include <math.h>
#include <iostream>
#include <string>

#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

#include "json.hpp"

using json = nlohmann::json;

static char help[] = "";

#define TYPEGEOM 1

int main(int ac, char **av) {

  std::ifstream i("file.json");
  json j;
  i >> j;

  std::cout << j;

#if defined (ENABLE_MPI)
  MPI_Init(0, NULL);

  // MPI information
  int nProcessors, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  int level = 4;
  double temperature = 0.8;
  double coef_penalization = 1.0e20;
  double tolerance_penalization = 0.2;
  Grid *g;
  LevelSet *ls;

  if (ac != 2) {
    std::cerr << av[0] << " 1 for ASphere, 2 for STL" << std::endl;
    exit(1);
  }

  int typeGeom = std::stoi(av[TYPEGEOM]);

  //Grid *g = new Grid(0.0, 0.0, 0.0, 1.0, 1.0/pow(2,level), 3);
  //ASphere *geo  = new ASphere(0.5, 0.5, 0.5, 0.2, 3);

  switch (typeGeom) {
  case 1: {
    g = new Grid(0.0, 0.0, 0.0, 1.0, 1.0/pow(2,level), 3);
    ASphere *geoSphere  = new ASphere(0.5, 0.5, 0.5, 0.2, 3);
    geoSphere->computeLevelSet(g);
    ls = new LevelSet(g);
    ls->addGeometry(geoSphere, "TAG");
    break;
  }
  default: {
    g = new Grid(-20, -10, -5, 40.0, 1, 3);
    STLGeometry *geoSTL = new STLGeometry("./Dragonite.stl", g);
    ls = new LevelSet(g);
    //geo->computeLevelSet(g);
    ls->addGeometry(geoSTL, "TAG");
  }
  }

  std::vector<double> phi0 = ls->getLevelSet("TAG");

  for (auto & element : phi0) {
    element = (element > tolerance_penalization) ? 0 : temperature;
  }

  PetscInitialize(&ac, &av, (char *)0, help);
  Laplacian *lap = (Laplacian *)LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, g);

  lap->setBCPositionCoef(BCPosition::OnInterface);
  lap->buildMatrix();
  lap->setDirichletCondition();
  lap->setRHS(phi0*=coef_penalization);
  lap->penalize(phi0/=temperature);

  bitpit::PiercedVector<double> sol = lap->solveLaplacian();

  std::vector<double> solv = PVtoV(sol, g);

  g->getVTK().setName("laplacian");

  g->getVTK().addData(
    "solution",
    bitpit::VTKFieldType::SCALAR,
    bitpit::VTKLocation::CELL,
    solv
    );

  g->getVTK().addData(
    "phi",
    bitpit::VTKFieldType::SCALAR,
    bitpit::VTKLocation::CELL,
    phi0
    );

  g->write();

  delete g;

#if defined (ENABLE_MPI)
  MPI_Finalize();
#endif
}
