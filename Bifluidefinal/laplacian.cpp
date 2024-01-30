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
#include <math.h>
#include <iostream>
#include <string>

#define TYPEGEOM 1

using namespace neos;

int main(int ac, char **av) {

  Neos_Init();

  int level = 2; // 4
  int dimension= 3;
  int maxr=3000000;
  double temperature = 0.8;
  double coef_penalization = 1.0e20;
  double tolerance_penalization = 0.1;
  Grid *g;
  LevelSet *ls;
  Geometry *geom;
  int geomId;

  if (ac != 2) {
    std::cerr << av[0] << " 1 for ASphere, 2 for STL" << std::endl;
    exit(1);
  }

  int typeGeom = std::stoi(av[TYPEGEOM]);

  //Grid *g = new Grid(0.0, 0.0, 0.0, 1.0, 1.0/pow(2,level), 3);
  //ASphere *geo  = new ASphere(0.5, 0.5, 0.5, 0.2, 3);

  switch (typeGeom) {
  case 1: {
    g = new Grid(0.0, 0.0, 0.0, 1.0, 1.0/pow(2,level), dimension);
    geom  = new ASphere(0.5, 0.5, 0.5, 0.2, dimension);
    //ACylinder *geoSphere  = new ACylinder(0.5, 0.5, 0.5, 0.2, dimension,2);
    break;
  }
  default: {
    //g = new Grid(-20, -10, -5, 20.0, 0.5, dimension);
    g = new Grid(-20, -20, -10, 50.0,   //40.0,
                 1, dimension);
    geom = new STLGeometry("./Dragonite.stl", g);
  }
  }

  ls = new LevelSet(g);
  geomId = ls->addGeometry(geom, "TAG");
  ls->compute(geomId);

  std::cerr << "Cells Count " << g->getCellCount() << std::endl;

  std::vector<double> phi0 = ls->getLevelSet("TAG");

  std::cerr << " raffinement " << g->getCellCount() << std::endl;

  PiercedVector<double> phi0PV = VtoPV(phi0, g);

  for (auto &cell : g->getCells()) {
    const long &id = cell.getId();

    if (0) {
      g->markCellForRefinement(id);
      //g->markCellForCoarsening(id);
    } else {
      if (g->getCell(id).isInterior()
          &&
          (phi0PV[id] < tolerance_penalization)
          &&
          (phi0PV[id] > 0.0)
          &&
          maxr-- > 0) {
        g->markCellForRefinement(id);
        //std::cerr << " raf " << maxr << std::endl;
      }
    }
  }

  std::cerr << " raffinement " << g->getCellCount() << std::endl;

  g->update();

  std::cerr << " raffinement " << g->getCellCount() << std::endl;

  ls->compute(geomId);
  phi0 = ls->getLevelSet("TAG");
  std::cerr << " phi  "  << std::endl;

  for (auto & element : phi0) {
    element = (element > tolerance_penalization) ? 0 : temperature;
  }

  Laplacian *lap = (Laplacian *)LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, g);

  lap->setBCPositionCoef(BCPosition::OnInterface);
  lap->buildMatrix();
  lap->setDirichletCondition();
  lap->setRHS(phi0*=coef_penalization);
  lap->penalize(phi0/=temperature);

  lap->solveLaplacian();

  PiercedVector<double> sol = lap->getSolutionPV();

  std::vector<double> solv = PVtoV(sol, g);

  g->setExportName("laplacian");

  g->addData(
    "solution",
    solv
    );

  g->addData(
    "phi",
    phi0
    );

  g->write();

  delete g;

  Neos_Finalize();

}
