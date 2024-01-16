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
#include "Neos.hpp"
#include "prediction-projection/Prediction.hpp"

#include <chrono>

//#include "patch_kernel.hpp"
#include "math.h"
#include <numeric>

#include "patch_kernel.hpp"
#include <cmath>


/**
 * This file contains an example of Poisson equation resolution.
 * The equation is:
 *                 \Delta u = f
 * with f = 4*cos(x+y)sin(x-y). Then u = cos(x+y)sin(x-y).
 * We Prescribe exact Dirichlet boundary conditions.
 **/


double u(NPoint pt, double t=0.)
{
  double x=pt[0], y=pt[1];
  return x*x +2*x;
}

double laplacianU(NPoint pt)
{
  double x=pt[0], y=pt[1];
  return 0;
}

using namespace neos;

int main(int ac, char **av) {

  Neos_Init(&ac,&av);

  /*----------------- Boundary conditions ------------------*/
  BoundaryConditions BC;
  BC.addCondition(0, "Dirichlet", u, Var::P);
  BC.addCondition(1, "Dirichlet", u, Var::P);
  BC.addCondition(2, "Dirichlet", u, Var::P);
  BC.addCondition(3, "Dirichlet", u, Var::P);


  auto start = std::chrono::high_resolution_clock::now();
  Grid *grid = new Grid(0, 0, 0, 1, 0.01, &BC, GRID_2D);

  /*----------------- Non uniform grid ------------------*/
  // UNCOMMENT ONE OF THE TWO FOR LOOP IF YOU WANT A NON UNIFORM GRID

  /*------------------- Refinment by boxes --------------*/
  //for (int k=0; k<3; k++)
  //{
  //	for (auto &cell : grid->getCells()) {
  //	    const long &id = cell.getId();
  //	    NPoint pt = grid->evalCellCentroid(id);
  //	    for (auto &v: pt)
  //	    {
  //		v = std::abs(v);
  //	    }
  //	    double xInf = *std::max_element(pt.begin(), pt.end());
  //
  //	    if( xInf <= (M_PI/4.+ 0.3 * k))
  //	    {
  //		grid->markCellForRefinement(id);
  //	    }
  //
  //	}
  //	grid->update(false,true);
  //}

  /*----------------------- Randomic refinment -----------------*/
  //for (int k=0; k<3; k++)
  //{
  //	for (auto &cell : grid->getCells())
  //	{
  //	    const long &id = cell.getId();
  //	    if (id%5 ==0)
  //		grid->markCellForRefinement(id);
  //
  //	}
  //	grid->update(false,true);
  //}


  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);

  std::cout<<"---------fin creation grille------"<<std::endl;
  std::cout << "Calculation time : " << duration.count() <<" seconds\n"<<std::endl;

  /*------------------ Build useful stencils--------------*/
  StencilBuilder stencils(grid);
  stencils.buildCellGradientStencil();
  stencils.buildInterfaceGradientStencil();
  
  


  /*---------------- Build FV Laplacian matrix ---------*/
  Laplacian *lap = LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, grid);
  lap->setStencils(&stencils);
  lap->toggleNeumann(false);
  int cellId,intId;
  PiercedVector<double> kappaCC, kappaFC, rhs;
  for (auto& cell: grid->getCells())
  { 
    cellId = cell.getId();
    NPoint octCenter = grid->evalCellCentroid(cellId);
    kappaCC.emplace(cellId);
    
      kappaCC[cellId] = 1./(1+octCenter[0]);
    
    
    rhs.emplace(cell.getId());
  }
  for (auto& inter: grid->getInterfaces())
  { 
    intId=inter.getId();
    NPoint octCenter = grid->evalInterfaceCentroid(intId);
    kappaFC.emplace(intId);
      
        kappaFC[intId] =1./(octCenter[0]+1); 

      
    
      
    
  }

  

  lap->buildFVMatrix(kappaCC,
                    kappaFC,
                    0,
                    Var::P);

  /*----------------------- Fill RHS ----------------------*/
  for (auto& cell: grid->getCells())
  {
    if (cell.isInterior())
    {
      int cellId = cell.getId();
      rhs[cellId] = -laplacianU(grid->evalCellCentroid(cellId));
    }
  }
  lap->addRHS(rhs);

  /* ----------------- Solve ---------------------- */
  lap->solveLaplacian();

  /*----------- Compare numerical and exact solution ----------*/
  std::vector<double> U = lap->getSolution();
  std::vector<double> Uex;
  double errInf(0.);
  double errL1(0.);

  for (auto &cell : grid->getCells())
  {
    const long &id = cell.getId();
    if (grid->getCell(id).isInterior())
    {
      Uex.push_back(u(grid->evalCellCentroid(id)));
    }
  }
  for (size_t i=0; i<Uex.size(); i++)
  {
    if ( std::abs(Uex[i] - U[i]) > errInf )
    {
      errInf = std::abs(Uex[i] - U[i]);
    }
    errL1 += std::abs(Uex[i] - U[i]);
  }

  /*---------------- Write VTK solutions ---------------*/
  grid->addData("U", U);

  grid->addData("Uex", Uex);

  grid->setExportName("LaplacianCos");
  grid->write();

  double minErrInf;
  double sumErrL1;
  double sumCells;
  double nCells = grid->nbCells();

  int rank;
  MPI_Reduce(&errInf,
             &minErrInf,
             1,
             MPI_DOUBLE,
             MPI_MAX,
             0,
             grid->getCommunicator());

  MPI_Reduce(&errL1,
             &sumErrL1,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             grid->getCommunicator());

  MPI_Reduce(&nCells,
             &sumCells,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             grid->getCommunicator());
  MPI_Comm_rank(grid->getCommunicator(), &rank);
  MPI_Barrier(grid->getCommunicator());


  double h = grid->getMinSize();

  if (rank == 0)
  {
    std::cout<<"Erreur: h:"<<h
             <<" L_{inf}: "<<minErrInf
             <<" L_1: "<<sumErrL1/sumCells<<std::endl;
  }

  delete grid;

  Neos_Finalize();

  return 0;
}
