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
#include "LaplacianFactory.hpp"
#include <math.h>
#include <iostream>
#include <string>

#define TYPEGEOM 1
#define MAXRAF 1

using namespace neos;

int main(int ac, char **av) {

  Neos_Init(&ac,&av);

  int level = 4; // 4
  int dimension= 3; // 3
  int maxr= std::stoi(av[MAXRAF]);
  double temperature_core = 0.8;
  double temperature_border = 0.2;
  double temperature_ambiant = 0.2;
  double coef_penalization = 1.0e20; // 1.0e20
  double tolerance_penalization = 0.1; // 0.1
  double coef_refinment = 0.2;
  double delta_t= 0.001; // 0.03; // time step
  int time_iteration= 90;
  int time_enable= 1; // 0: disable, 1: enable time step
  int border_enable= 0; // 0 disable, 1: enable border penalisation

  double h=1;  // grid step -

  Grid *g;
  LevelSet *ls;

  if (!time_enable)
  {
    time_iteration= 1;
  }

  if (ac != 2) {
    std::cerr << av[0] << " 1 for ASphere, 2 for STL" << std::endl;
    exit(1);
  }

  int typeGeom = std::stoi(av[TYPEGEOM]);

  int geoId = -1;

  switch (typeGeom) {
    case 1: {
      g = new Grid(0.0, 0.0, 0.0, 1.0, 1.0/pow(2,level), dimension);
      g->update(false,true);
      ASphere *geoSphere  = new ASphere(0.5, 0.5, 0.5, 0.2, dimension);
      geoSphere->computeLevelSet(g);
      ls = new LevelSet(g);
      geoId = ls->addGeometry(geoSphere, "TAG");
      break;
    }
    default: {
      g = new Grid(-20, -20, -10, 50.0, 2.0, dimension);
      STLGeometry *geoSTL = new STLGeometry();
      geoSTL->addSTL("./Dragonite.stl", g);
      ls = new LevelSet(g);
      geoSTL->computeLevelSet(g);
      geoId = ls->addGeometry(geoSTL, "TAG");
      delta_t= 0.1;
      coef_refinment = 20.0;
    }
  }

  std::vector<double> phi0 = ls->getLevelSet("TAG");

  int nbRef = g->refineFromRangeVector(phi0,0,coef_refinment);
  std::cout << nbRef << " Cells Refined" << '\n';
  Geometry *geom = ls->getGeometry(geoId);
  geom->computeLevelSet(g);
  phi0 = ls->getLevelSet("TAG");

  double h_h= pow(h,dimension);
  double h_h_sur_delta_t= h_h/delta_t;

  std::cerr << " h : " << h << " h h " << h_h << " delta_t " << delta_t << " h_h_sur_delta_t " << h_h_sur_delta_t << std::endl;

  std::vector<double> pen = phi0;
  std::vector<double> initial_temperature_RHS = phi0;

  for (auto &element : pen) {
    if (element > tolerance_penalization) {
      element = 0.0; // no penalization
    } else {
      element = 1.0; // penalization inside
    }
  }
  std::vector<double> fixed_temperature_RHS = pen * temperature_core;

  for (auto &element : initial_temperature_RHS) {
    if (element > tolerance_penalization) {
      element = temperature_ambiant;
    } else {
      element = temperature_core;
    }
  }

  std::vector<double> sol_v= initial_temperature_RHS;

  if (0) {
    for (auto &element : sol_v) {
      element = 0.0;
    }
  }

  std::vector<std::array<double, 3> > gradvtk;

  g->setExportName("time-laplacian");
  g->setCounter();

  Laplacian *lap = (Laplacian *)LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, g);

  lap->setBCPositionCoef(BCPosition::OnInterface);
  lap->buildMatrix();

  lap->penalize(pen * coef_penalization * (-1.0) * h_h - time_enable * h_h_sur_delta_t);

  for (int i=0; i<time_iteration; i++) {
    lap->setRHS( fixed_temperature_RHS * coef_penalization * (-1.0) * h_h - time_enable * h_h_sur_delta_t * sol_v); // RH

    double norm = lap->solveLaplacian();
    sol_v = lap->getSolution();
    std::cout << "Iter " << i << " Norm " << norm << std::endl;

    g->addData("solution",sol_v);
    g->addData("phi",phi0);

    g->write();
  }
  delete g;

  Neos_Finalize();
}
