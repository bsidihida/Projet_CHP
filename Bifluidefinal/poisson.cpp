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
#include "InterpolatorFactory.hpp"

#include "patch_kernel.hpp"

#define TYPEGEOM 1
#define MAXRAF 1

using namespace neos;

int main(int ac, char **av) {

  Neos_Init(&ac, &av);

  #if defined (ENABLE_MPI)
  // MPI information
  int nProcessors, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  int level = 3;
  int dimension= 2;
  int maxr= std::stoi(av[MAXRAF]);
  Grid *g;

  if (ac != 2) {
    std::cerr << av[0] << " nb raffinement missing" << std::endl;
    exit(1);
  }

  g = new Grid(-1.0, -1.0, 0.0, 2.0, 2.0/pow(2,level), dimension);

  std::cerr << " raffinement " << g->getCellCount() << std::endl;

  // raffinement
  for (auto &cell : g->getCells())
  {
    const long &id = cell.getId();
    std::array<double, 3> ic;

    if (g->getCell(id).isInterior()) {
      ic = g->evalCellCentroid(id);
      double x=ic[0];
      double y=ic[1];


      double val = -2*pow(3.141592653589,2)*sin(3.141592653589*x)*sin(3.141592653589*y);

      if ( val >0.5 && maxr>=0 ) {
        g->markCellForRefinement(id);
      }
    }
  }

  g->update(false,true);
  std::cerr << " raffinement count: " << g->getCellCount() << " nb: " << g->nbCells() << std::endl;

  std::vector<double> vg0(g->nbCells());
  std::vector<double> vg1(g->nbCells());
  std::vector<double> afafb(g->nbCells());
  std::vector<double> xy(g->nbCells());
  std::vector<double> refid(g->nbCells());
  std::vector<double> sol;
  std::vector<std::array<double, 3> > gradvtk;

  double lpx=0.1,lpy=0.1,lpz=0.1;
  NPoint nlp={0.3,0.3,0.5};

  std::cerr << "celule id: " << g->locatePoint(nlp) << std::endl;
  std::cerr << "celule id: " << g->locatePoint({lpx,lpy,lpz})  << std::endl;

  long cid= g->locatePoint(nlp);

  std::vector<long> cidnei;
  cidnei= g->findCellNeighs(cid);

  std::vector<std::array<double, 3> > xref;
  std::vector<double> val_ref;

  double valnawouac= 0.5;

  for (auto &nei : cidnei)
  {
    std::cerr << "nei: " << nei << std::endl;
    xref.push_back(g->evalCellCentroid(nei));
    val_ref.push_back(valnawouac++);
  }

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::DISTANCEWEIGHTED);

  double res = interpo->computeInterpolation(nlp,
                                             xref,
                                             val_ref);

  std::cerr << "interpol= " << res << std::endl;


  long count = 0;
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

      g0 = -2*pow(3.141592653589,2)*sin(3.141592653589*x)*sin(3.141592653589*y);
      g1 = sin(3.141592653589*x)*sin(3.141592653589*y);

      vg0[count]= g0;
      vg1[count]= g1;

      xy[count]= int(x*1000.0)+(int(y*1000.0)/1000.0); // xy = x,y

      afafb[count]= g0;

      refid[count]= id;

      count++;
    }
  }

  Laplacian *lap = (Laplacian *)LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, g);

  lap->setBCPositionCoef(BCPosition::OnInterface);
  lap->buildMatrix();
  lap->setDirichletCondition();
  lap->setRHS(afafb);
  lap->solveLaplacian();
  sol = lap->getSolution();

  MPI_Barrier(g->getCommunicator());


  std::cerr << "XXX getGrad..." << std::endl;
  gradvtk= lap->getGradientFromSolution();

  std::vector<double> gradvtk_x(gradvtk.size());
  std::vector<double>::iterator iG = gradvtk_x.begin();

  for (auto &x : gradvtk)
  {
    *iG= x[0];
    ++iG;
  }


  std::vector<double> gradlap_x(g->nbCells());
  std::vector<double> gradlap_y(g->nbCells());

  std::vector<double> gradvtklolox(g->nbCells());

  count=0;
  for (auto &cell : g->getCells())
  {
    const long &id = cell.getId();
    std::array<double, 3> ic;

    if (g->getCell(id).isInterior()) {
      ic = g->evalCellCentroid(id);

      gradvtklolox[count]= gradvtk[count][0];

      gradlap_x[count]= gradvtk[count][0];
      gradlap_y[count]= gradvtk[count][1];


      count++;
    }
  }
#ifdef NONO
  Gradient hessien(g->getDimension(), g);
  std::vector<std::array<double, 3> > hess_x= hessien.computeLSGradient(gradvtk_x);

  std::vector<std::array<double, 3> > hess_gradvtklolox= hessien.computeLSGradient(gradvtklolox);
  std::vector<std::array<double, 3> > hesslap_x= hessien.computeLSGradient(vg0);
  std::vector<std::array<double, 3> > hesslap_y= hessien.computeLSGradient(gradlap_y);
#endif
  g->setExportName("poisson-afaf");
#ifdef NONO
  g->addData("solution",sol);

  g->addData(
    "hessien_x"
    hess_x
    );

#endif
  g->addData("GradLapx",gradlap_x);
  g->addData("add_gradvtklolox",gradvtklolox);

#ifdef NONO
  g->addData(
    "LapHessien_x",
    hesslap_x
    );

  g->addData(
    "LapHessien_y",
    hesslap_y
    );


  g->addData(
    "gradientVTK",
    gradvtk
    );

  g->addData(
    "hess_gradvtklolox",
    hess_gradvtklolox
    );
#endif

  g->write();

  delete g;

  Neos_Finalize();
}
