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
#include <chrono>

#include "Neos.hpp"
#include "Transport.hpp"
#include "CommonTools.hpp"
#include "SimulationParameters.hpp"
#include "Elasticity.hpp"
#include "prediction-projection/Projection.hpp"
#include "prediction-projection/Prediction.hpp"

double F(double t)
{
  double nu = 0.01;
  return exp(-2*nu*t);
}

double taylorGreen_Ux(NPoint pt,double t=0.)
{
  return sin(pt[0]) * cos(pt[1]) * F(t);
}

double taylorGreen_Uy(NPoint pt,double t=0.)
{
  return -cos(pt[0]) * sin(pt[1]) * F(t);
}


double phi0(NPoint pt,double yc)
{
  
  return  (pt[1]-yc);
  
  
}
double bord(NPoint pt,double t=0.)
{
  return 0;
}
double un(NPoint pt,double t=0){

  return 1;
}





using namespace neos;

/**
 * TaylorGreen.cpp
 * Example to show how to use neos to compute a taylor green vortex in 2D
 * We have comparison with analytical solution at the end for Tmax = 1.
 **/

int main(int ac, char **av)
{
  Neos_Init(&ac,&av);

  /*------------------------ Load Parameters----------------------*/
  // FIXME: CHANGE FILE NAME HERE
  SimulationParameters param;
  param.readFromInputFile("InputTG.dat");

  auto start = std::chrono::high_resolution_clock::now();

  Solution sol, solPrev, solNext;

  /*----------------- Boundary conditions ------------------*/
  BoundaryConditions BC;
  BC.addCondition(0, "Dirichlet", bord, Var::Ux);
  BC.addCondition(1, "Dirichlet", bord, Var::Ux);
  BC.addCondition(2, "Dirichlet", bord, Var::Ux);
  BC.addCondition(3, "Dirichlet", bord, Var::Ux);

  BC.addCondition(0, "Dirichlet", bord, Var::Uy);
  BC.addCondition(1, "Dirichlet", bord, Var::Uy);
  BC.addCondition(2, "Dirichlet", bord, Var::Uy);
  BC.addCondition(3, "Dirichlet", bord, Var::Uy);

  BC.addCondition(0, "Neumann", taylorGreen_Uy, Var::P);
  BC.addCondition(1, "Neumann", taylorGreen_Uy, Var::P);
  BC.addCondition(2, "Neumann", taylorGreen_Uy, Var::P);
  BC.addCondition(3, "Dirichlet", taylorGreen_Uy, Var::P);

  Grid     *grid = new Grid(0,0,0, 2*M_PI, 0.1, &BC, GRID_2D);
  
  PiercedVector<double> phiCC,phiFC;
  PiercedVector<double> divg, divg2;
  PiercedVector<double> phiTBuff;
  /*----------------- Compute initial velocity --------------*/
  for (auto &cell: grid->getCells())
  {
    int cellId = cell.getId();

    sol.Ux.emplace(cellId);
    sol.Uy.emplace(cellId);
    phiCC.emplace(cellId);
    phiTBuff.emplace(cellId);

    phiCC[cellId] = phi0(grid->evalCellCentroid(cellId),3.);

    sol.Ux[cellId] = un(grid->evalCellCentroid(cellId));
    sol.Uy[cellId] = un(grid->evalCellCentroid(cellId));
    
    sol.VelocityCC.emplace(cellId);
    sol.VelocityCC[cellId] = {sol.Ux[cellId],sol.Uy[cellId],0.};

    sol.FctInd.emplace(cellId);
    sol.FctInd[cellId]     = 1000.;

    solNext.FctInd.emplace(cellId);
    solNext.FctInd[cellId]     = 1000.;

    sol.Pressure.emplace(cellId);
    sol.Pressure[cellId] = 0.;

    solNext.Pressure.emplace(cellId);
    solNext.Pressure[cellId] = 0.;

    solNext.Ux.emplace(cellId);
    solNext.Uy.emplace(cellId);
    solNext.VelocityCC.emplace(cellId);
  }


  for (auto& inter:grid->getInterfaces())
  {

  phiFC.emplace(inter.getId());
  phiFC[inter.getId()] = phi0(grid->evalInterfaceCentroid(inter.getId()),3.);
    
  }

  /*--------------- Write solution at time 0 ---------- */
  auto velCC = PVtoV(sol.VelocityCC, grid);
  auto pressure = PVtoV(sol.Pressure, grid);
  grid->setExportName("TaylorGreen_0");
  grid->addData("Vel", velCC);
  grid->addData("pressure", pressure);
  grid->write();

  /*----------------- Projection at face center ----------*/
  StencilBuilder stencils(grid);
  stencils.buildCellGradientStencil();
  stencils.buildInterfaceGradientStencil();
  stencils.buildCCToFCStencils();

  Projection proj(grid, &stencils, &param);

  sol.VelocityFC = proj.FaceCenterProjection(sol.VelocityCC);
  solNext.VelocityFC = proj.FaceCenterProjection(sol.VelocityCC);

  /*---------------- Solution initialization ------------ */
  solPrev.VelocityFC = sol.VelocityFC;
  solPrev.VelocityCC = sol.VelocityCC;


  solPrev.Ux = sol.Ux;
  solPrev.Uy = sol.Uy;

  solPrev.Pressure = sol.Pressure;

  solPrev.t = 0.;
  sol.t = 0;
  solNext.t = 0;

  solPrev.dt = 0.005;
  sol.dt = 0.005;
  solNext.dt = 0.005;

  solPrev.FctInd = sol.FctInd;

  int iter(0);
  Prediction prediction(grid, &stencils, &param);
  Transport trspt(grid, &stencils);

  /*-------------------  Time loop -------------------*/
  while (solNext.t < 0.1)


  {

  double t=0.;
   double dt = 0.0001;
   while (t<solNext.dt)
   {

  //   // Compute phi at time n+1 using RK2 scheme
    divg = trspt.computeWithSecondOrderFV(sol.VelocityFC,
                                          phiCC,
                                          t,
                                          Var::P);
   


    

  //   // Phi^{n+1./2}
  
  
    t += dt/2.;
    for (auto &id: divg.getIds())
     {
       phiTBuff[id] = phiCC[id] - 0.5*dt*divg[id];
     }


    divg2 = trspt.computeWithSecondOrderFV(sol.VelocityFC,
                                          phiTBuff,
                                          t,
                                          Var::P);
  

     t += dt/2.;
    // Phi^{n+1}
    for (auto &id: divg2.getIds())
    {
      phiCC[id] -= dt * divg2[id];
    }

   }

   

    for (auto& inter:grid->getInterfaces())
  {
    

    const long& interId = inter.getId();
    const auto& owners  = inter.getOwnerNeigh();

    phiFC[interId]=phiCC[owners[0]];
    printf("phiFCest%f\n",phiFC[interId]);
   
   

    
  }


    


 







    

    prediction.ComputePredictionStep(solPrev,
                                     sol,
                                     solNext,phiFC,phiCC);
    printf("test\n");

    proj.ComputeProjectionStep(solNext,
                               sol,phiFC,phiCC);

    solNext.t += 0.005;

    CommonTools::update(solPrev, sol, solNext);

    iter += 1;
    if (iter%10 == 0)
    {
      velCC = PVtoV(sol.VelocityCC, grid);
      pressure = PVtoV(sol.Pressure, grid);
      grid->setExportName("TaylorGreen_" + std::to_string(iter));
      grid->addData("Vel", velCC);
      grid->addData("pressure", pressure);
      grid->write();
    }

    std::cout<<solNext.t<<std::endl;
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << "Calculation time : " << duration.count() <<"seconds\n"<<std::endl;

  /* --------------------- Error computation ---------------------- */
  double errInfX(0.), errInfY(0.);
  double errL1X(0.), errL1Y(0.);

  for (auto &cell: grid->getCells())
  {
    if (cell.isInterior())
    {
      const long& cellId = cell.getId();
      NPoint pt = grid->evalCellCentroid(cellId);

      std::array<double,3> velocityEx =
      {taylorGreen_Ux(pt, 0.1), taylorGreen_Uy(pt, 0.1), 0};

      if (std::abs(solNext.VelocityCC[cellId][0] - velocityEx[0])
          > errInfX)
        errInfX = std::abs(solNext.VelocityCC[cellId][0] - velocityEx[0]);

      if (std::abs(solNext.VelocityCC[cellId][1] - velocityEx[1])
          > errInfY)
        errInfY = std::abs(solNext.VelocityCC[cellId][1] - velocityEx[1]);

      errL1X += std::abs(solNext.VelocityCC[cellId][0] -
                         velocityEx[0]);
      errL1Y += std::abs(solNext.VelocityCC[cellId][1] -
                         velocityEx[1]);

    }
  }

#if defined (ENABLE_MPI)
  double minErrInfX, sumErrL1X;
  double minErrInfY, sumErrL1Y;
  double sumCells;
  double nCells = grid->nbCells();
  int rank;
  MPI_Reduce(&errInfX,
             &minErrInfX,
             1,
             MPI_DOUBLE,
             MPI_MAX,
             0,
             grid->getCommunicator());

  MPI_Reduce(&errInfY,
             &minErrInfY,
             1,
             MPI_DOUBLE,
             MPI_MAX,
             0,
             grid->getCommunicator());

  MPI_Reduce(&errL1X,
             &sumErrL1X,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             grid->getCommunicator());

  MPI_Reduce(&errL1Y,
             &sumErrL1Y,
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
    std::cout<<"Erreur sur X: h:"<<h
             <<" L_{inf}: "<<minErrInfX
             <<" L_1: "<<sumErrL1X/sumCells<<std::endl;
    std::cout<<"Erreur sur Y: h: "<<h
             <<" L_{inf}: "<<minErrInfY
             <<" L_1: "<<sumErrL1Y/sumCells<<std::endl;
  }

#endif
  delete grid;

  Neos_Finalize();

  std::cout<<"called petsc finalize"<<std::endl;
  return 0;
}
