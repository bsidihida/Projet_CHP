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
//#include "bitpit.hpp"
#include "prediction-projection/Projection.hpp"
#include <chrono>

#include "patch_kernel.hpp"
#include <cmath>

using namespace neos;

double phi(NPoint pt, double t=0.)
{
  double x = pt[0], y=pt[1];
  double R = sqrt(x*x + y*y);
  double f = tanh(R)/(cosh(R) * cosh(R));
  double temp = (f *t) / (0.385 * R);
  return -tanh(0.5 * y *cos(temp) - 0.5 * x * sin(temp));
}

double phi2(NPoint pt, double t=0.)
{
  double x = pt[0], y=pt[1];
  double R = sqrt(x*x + y*y);
  double f = tanh(R)/(cosh(R) * cosh(R));
  double temp = f / (0.385 * R);
  return temp * -y;
}

double phi3(NPoint pt, double t=0.)
{
  double x = pt[0], y=pt[1];
  double R = sqrt(x*x + y*y);
  double f = tanh(R)/(cosh(R) * cosh(R));
  double temp = f / (0.385 * R);
  return temp * x;
}

NPoint T(double x, double y, int octId, Grid *g)
{
  NPoint octCenter = g->evalCellCentroid(octId);
  double h = g->evalCellSize(octId);

  return {h/2. * x + octCenter[0], h/2. * y + octCenter[1], 0};
}

double computeIntegralPhi(int octId, Grid *g, double t=0.)
{
  double h = g->evalCellSize(octId);

  return h*h/4.
         * ( phi(T(sqrt(1./3), sqrt(1./3), octId, g),t)
             +   phi(T(-sqrt(1./3), sqrt(1./3), octId, g),t)
             +   phi(T(-sqrt(1./3), -sqrt(1./3), octId, g),t)
             +   phi(T(sqrt(1./3), -sqrt(1./3), octId, g),t ));

}

int main(int ac, char **av) {

  Neos_Init(&ac,&av);


  /*------------- Add Boundary conditions ----------*/
  BoundaryConditions BC;

  // BC for phi function
  BC.addCondition(0, "Dirichlet", phi, Var::P);
  BC.addCondition(1, "Dirichlet", phi, Var::P);
  BC.addCondition(2, "Dirichlet", phi, Var::P);
  BC.addCondition(3, "Dirichlet", phi, Var::P);

  // BC for velocity
  BC.addCondition(0, "Dirichlet", phi2, Var::Ux);
  BC.addCondition(1, "Dirichlet", phi2, Var::Ux);
  BC.addCondition(2, "Dirichlet", phi3, Var::Uy);
  BC.addCondition(3, "Dirichlet", phi3, Var::Uy);

  // Init grid
  Grid *g = new Grid(-1, -1, 0, 2, 0.05, &BC, GRID_2D);
  g->update(false,true);
  std::cout<<"---------fin creation grille------"<<std::endl;

  // Define velocity at cell center, phi function, vel at face
  // center and divergence for finite volume transport
  PiercedVector<bitpit::darray3> velCC;
  PiercedVector<double> phiT;
  PiercedVector<double> phiTBuff;
  PiercedVector<double> velFC;
  PiercedVector<double> divg, divg2;
  StencilBuilder stencils(g);
  stencils.buildCellGradientStencil();
  stencils.buildInterfaceGradientStencil();

  Projection proj(g);
  Transport trspt(g, &stencils);
  //Transport trspt(g);
  int cellId;

  /*----------- Init velocity and phi function ------------*/
  for (auto &cell: g->getCells())
  {
    cellId = cell.getId();
    velCC.emplace(cellId);

    double cellVolume = g->evalCellVolume(cellId);
    NPoint octCenter = g->evalCellCentroid(cellId);
    double R = sqrt(octCenter[0]*octCenter[0]
                    + octCenter[1]*octCenter[1]);
    double f = tanh(R)/(cosh(R) * cosh(R));
    double temp = f / (0.385 * R);

    velCC[cellId] = {-octCenter[1], octCenter[0], 0};
    velCC[cellId] *= temp;

    phiT.emplace(cellId);
    phiTBuff.emplace(cellId);
    phiT[cellId]  = 1./cellVolume * computeIntegralPhi(cellId, g);;
  }

  /*--------- Projection of velocity at face center ----------*/
  velFC = proj.FaceCenterProjection(velCC);

  double t=0.;
  double dt = 0.001;
  auto i = 0u;

  /*--------------------- Time loop -------------------*/
  auto start = std::chrono::high_resolution_clock::now();
  g->setExportName("transportAtmo"+std::to_string(i));
  g->setCounter();
  while (t<0.5)
  {

    // Compute phi at time n+1 using RK2 scheme
    divg = trspt.computeWithSecondOrderFV(velFC,
                                          phiT,
                                          t,
                                          Var::P);

    // Write outputs
    std::vector<double> divgV = PVtoV(divg, g);
    std::vector<double> phiTV = PVtoV(phiT, g);
    if (i%10 == 0)
    {
      g->addData("phi", phiTV);
      g->addData("div", divgV);
      g->write();
    }

    // Phi^{n+1./2}
    t += dt/2.;
    for (auto &id: divg.getIds())
    {
      phiTBuff[id] = phiT[id] - 0.5*dt*divg[id];
    }


    divg2 = trspt.computeWithSecondOrderFV(velFC,
                                          phiTBuff,
                                          t,
                                          Var::P);

    t += dt/2.;
    // Phi^{n+1}
    for (auto &id: divg2.getIds())
    {
      phiT[id] -= dt * divg2[id];
    }

    i++;
    std::cout<<"\rItération: "<<i<<std::flush;
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << "Calculation time : " << duration.count()
            <<" seconds\n"<<std::endl;

  delete g;

  Neos_Finalize();

  return 0;
}
