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
#define MAXRAF 1

using namespace neos;

int main(int ac, char **av) {

  Neos_Init(&ac,&av);

  int level = 6; // 4
  int dimension= 2; // 3
  int time_iteration= 20;
  int maxr= std::stoi(av[MAXRAF]);

  Grid *g;

  g = new Grid(0.0, 0.0, 0.0, 1.0, 1.0/pow(2,level), dimension);

  // raffinement
  for (auto &cell : g->getCells())
  {
    const long &id = cell.getId();
    std::array<double, 3> ic;

    if (g->getCell(id).isInterior()) {
      ic = g->evalCellCentroid(id);
      double x=ic[0];
      double y=ic[1];
      if (x*x<y && maxr-->0) {
        g->markCellForRefinement(id);
        //g->markCellForCoarsening(id);
      }
    }
  }

  g->update(false,true);
  std::cerr << " raffinement count: " << g->getCellCount() << " nb: " << g->nbCells() << std::endl;

  std::vector<double> sol;
  std::vector<std::array<double, 3> > gradvtk;
  std::vector<double> afafb(g->nbCells());
  std::vector<double> vg0(g->nbCells());
  std::vector<double> vg1(g->nbCells());

  g->setExportName("time-poisson");
  g->setCounter();

  Laplacian *lap = (Laplacian *)LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, g);

  lap->setBCPositionCoef(BCPosition::OnInterface);
  lap->buildMatrix();
  lap->setDirichletCondition();


  for (int i=0; i<time_iteration; i++) {

    int count= 0;
    for (auto &cell : g->getCells())
    {
      const long &id = cell.getId();
      std::array<double, 3> ic;

      if (g->getCell(id).isInterior()) {

        ic = g->evalCellCentroid(id);

        double x=ic[0];
        double y=ic[1];
        double g0=0;
        double g1=0;

        g0= exp(-10*( pow(x+(i/(time_iteration*2.0))-0.6,2)+pow(y-(i/(time_iteration*2.0))-0.6,2) )) + 0.001;
        g1= exp(-10*( pow(x-0.3,2)+pow(y-0.3,2) )) + 0.001;

        afafb[count]= (g0-g1)/g1;
        vg0[count]= g0;
        vg1[count]= g1;

        count++;
      }
    }
    lap->setRHS( afafb );
    lap->solveLaplacian();
    sol = lap->getSolution();
    gradvtk= lap->getGradientFromSolution();

    // calcul du Hessian
    std::vector<double> gradvtk_x(gradvtk.size());
    std::vector<double>::iterator iG = gradvtk_x.begin();
    for (auto &x : gradvtk)
    {
      *iG= x[0];
      ++iG;
    }

    Gradient hessien(g->getDimension(), g);
    std::vector<std::array<double, 3> > hess_x= hessien.computeLSGradient(gradvtk_x);

    g->addData("solution",sol);
    g->addData("afafb",afafb);
    g->addData("g0",vg0);
    g->addData("g1",vg1);
    g->addData("gradientVTK", gradvtk);
    g->addData("hessien_x",hess_x);
    g->write();
  }
  delete g;

  MPI_Finalize();

}
