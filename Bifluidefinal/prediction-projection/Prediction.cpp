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

/**
 * @file   Prediction.cpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Mon Sep  9 13:42:40 2019
 *
 * @brief  This file contains Prediction class
 *
 * @copyright Inria
 */

#include "Prediction.hpp"
#include "LaplacianFactory.hpp"
#include "Gradient.hpp"
#include "Transport.hpp"
#include <chrono>
#include <numeric>
#include <cmath>


double He3(double phi, double epsilon) {
    if (phi < -epsilon) {
        return 0.0;
    } else if (std::abs(phi) <= epsilon) {
        return 0.5 * (1 + phi / epsilon + sin(M_PI * phi / epsilon) / M_PI);
    } else {
        return 1.0;
    }
}
double rho3(double phi,double epsilon,double rho1,double rho2)
{
  return rho1+(rho2-rho1)*He3(phi,epsilon);
}



double mu3(double phi,double epsilon,double mu1,double mu2)
{
  return mu1+(mu2-mu1)*He3(phi,epsilon);
  
}


// double rho(NPoint pt,double rho1,double rho2)
// {
//   if(pt[1]>3)
//   return rho1;
//   else{

//     return rho2;
//   }
// }

// double mu(NPoint pt,double mu1,double mu2)
// {
//   if(pt[1]>3)
//   return mu1;
//   else{

//     return mu2;
//   }
// }
double rho1=1000.,rho2=1.,mu1=0.001,mu2=0.01,epsilon= 0;
namespace neos {

UserDataComm<double>* Prediction::_pvComm = NULL;
int Prediction::_pvComm_only_one_init_dirty=0;

Prediction::Prediction(Grid* grid,
                       StencilBuilder* stencils,
                       SimulationParameters* param ) : _grid(grid),

  _param(param),
  _hasParam(param != NULL),
  _stencils(stencils),
  _rhocc({}),
  _mufc({})
{
  lapU = LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, grid);
  lapV = LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, grid);
  _mapGlobal = _grid->getCellGlobalMap();



  for (auto& cell: _grid->getCells())
  {
    _rhocc.emplace(cell.getId());
   
  }
 for (auto& inter:_grid->getInterfaces())
  {

    _mufc.emplace(inter.getId());
   
  }

  if (!_hasParam)
  {
    _param = new SimulationParameters();
  }
}

Prediction::~Prediction()
{
  if(!_hasParam)
  {
    delete _param;
  }
}

void Prediction::ComputePredictionStep(Solution& solPrev,
                                       Solution& sol,
                                       Solution& solNext,PiercedVector<double> phiFC, PiercedVector<double> phiCC)
{
  for (auto& cell: _grid->getCells())
  {
  
    _rhocc[cell.getId()]=rho3(phiCC[cell.getId()], epsilon, rho1, rho2);
  }
  for (auto& inter:_grid->getInterfaces())
  {

  
    _mufc[inter.getId()]=mu3(phiFC[inter.getId()], epsilon, mu1, mu2);
  }


  solNext.VelocityCC.fill({0., 0., 0.});
  solNext.Ux.fill(0.);
  solNext.Uy.fill(0.);
  // Coeff for Gear scheme
  // FIXME: Add a class to handle time scheme.
  double c_n(0.), c_nm1(0.), c_f(0.);

  c_n = 1. + solNext.dt * solNext.dt / (sol.dt * (sol.dt + 2 * solNext.dt));
  c_nm1 = -solNext.dt * solNext.dt / (sol.dt * (sol.dt + 2*solNext.dt));
  c_f = solNext.dt * (solNext.dt + sol.dt) / (sol.dt + 2*solNext.dt);


  // Coeff for Adams-Bashforth scheme
  double cAB_n(0.), cAB_nm1(0.);

  cAB_n = 1. + solNext.dt/sol.dt;
  cAB_nm1 = -solNext.dt/sol.dt;

  // FIXME: it may be better to define an augmented matrix ?
  // Here we have a laplacian for X and Y composant.
  // FIXME: Only for 2D problems

  lapU->setStencils(_stencils);
  lapV->setStencils(_stencils);
  lapU->toggleNeumann(false);
  lapV->toggleNeumann(false);
  Transport trspt(_grid, _stencils);

  PiercedVector<double> kappaCC, kappaFC;
  for (auto& cell: _grid->getCells())
  {
    kappaCC.emplace(cell.getId());
    kappaCC[cell.getId()] = -0;//_param->getFluidKinematicViscosity() * c_f;//0;
  }
  for (auto& inter: _grid->getInterfaces())
  {
    kappaFC.emplace(inter.getId());
    kappaFC[inter.getId()] = 1.;
  }

  // FIXME: Separate build of matrix and vector.
  // The matrix could be computed once at each mesh refinment instead of at
  // each iteration (if there is no penalization).
  lapU->zeroMatrix();
  lapU->zeroRHS();
  lapU->buildFVMatrix(kappaCC,
                     kappaFC,
                     solNext.t + solNext.dt,
                     Var::Ux);
  lapV->zeroMatrix();
  lapV->zeroRHS();
  lapV->buildFVMatrix(kappaCC,
                     kappaFC,
                     solNext.t+solNext.dt,
                     Var::Uy);


  double eps = 1e-12;
  double coeffPen = c_f/eps;

  // Penalize for Ux and Uy
  lapU->penalizeAtOrder1(solNext.FctInd,
                        coeffPen);
  lapV->penalizeAtOrder1(solNext.FctInd,
                        coeffPen);

  // Compute pressure gradient
  Gradient grad(_grid->getDimension(), _grid, _stencils);
  PiercedVector<std::array<double, 3> > gradPressure
    = grad.computeFVGradient(sol.Pressure,
                             solNext.t,
                             Var::P);
  










const auto gradUx = grad.computeFCLSGradient(sol.Ux,solNext.t,Var::Ux);
const auto gradUy = grad.computeFCLSGradient(sol.Uy,solNext.t,Var::Uy);


PiercedVector<std::array<double, 2> > divergence;
double ier=0;
//ici j'ai besoin des gradient aux interfaces pas juste du gradient

  for (auto& inter:_grid->getInterfaces())
  {
    

    const long& interId = inter.getId();
    const auto& owners  = inter.getOwnerNeigh();
    const auto& normal = _grid->evalInterfaceNormal(interId);
    bool isVertical = std::abs(normal[0]) > std::abs(normal[1]);
   
    double gradientUx =  gradUx[interId] ;
    double gradientUy = gradUy[interId];

    if (_grid->getCell(owners[0]).isInterior())
    {
      if (!divergence.exists(owners[0]))
        divergence.emplace(owners[0]);
     

      divergence[owners[0]][0] += _mufc[inter.getId()]*(gradientUx)
                               
                               * _grid->evalInterfaceArea(interId)
                               / _grid->evalCellVolume(owners[0]);
      divergence[owners[0]][1] += _mufc[inter.getId()]*(gradientUy)
                               
                               * _grid->evalInterfaceArea(interId)
                               / _grid->evalCellVolume(owners[0]);
   

    

    if (owners[1] != bitpit::Element::NULL_ID)
    {
      if (_grid->getCell(owners[1]).isInterior())
      {
        if (!divergence.exists(owners[1]))
          divergence.emplace(owners[1]);

        divergence[owners[1]][0] -= _mufc[inter.getId()]*(gradientUx)
                                
                                 * _grid->evalInterfaceArea(interId)
                                 / _grid->evalCellVolume(owners[1]);
        divergence[owners[1]][1] -= _mufc[inter.getId()]*(gradientUy)
                                 
                                 * _grid->evalInterfaceArea(interId)
                                 / _grid->evalCellVolume(owners[1]);
   

      }
    }
    }
    
    


    
  }
                        

  // Transport all values that needed to be.
  // Ux^{n}-Uy^{n}-Ux{n-1}-Uy^{n-1}
  std::vector<PiercedVector<double> > valToTransport;
  std::vector<PiercedVector<double> > velsFC;
  std::vector<Var> types;
  std::vector<double> times;


  valToTransport.push_back(sol.Ux);
  velsFC.push_back(sol.VelocityFC);
  types.push_back(Var::Ux);
  times.push_back(sol.t);

  valToTransport.push_back(sol.Uy);
  velsFC.push_back(sol.VelocityFC);
  types.push_back(Var::Uy);
  times.push_back(sol.t);

  if (_grid->getDimension()>2)
  {
    valToTransport.push_back(sol.Uz);
    velsFC.push_back(sol.VelocityFC);
    types.push_back(Var::Uz);
    times.push_back(sol.t);
  }

  valToTransport.push_back(solPrev.Ux);
  velsFC.push_back(solPrev.VelocityFC);
  types.push_back(Var::Ux);
  times.push_back(solPrev.t);

  valToTransport.push_back(solPrev.Uy);
  velsFC.push_back(solPrev.VelocityFC);
  types.push_back(Var::Uy);
  times.push_back(solPrev.t);
  if (_grid->getDimension()>2)
  {
    valToTransport.push_back(solPrev.Uz);
    velsFC.push_back(solPrev.VelocityFC);
    types.push_back(Var::Uz);
    times.push_back(solPrev.t);
  }

  std::vector<PiercedVector<double> > trsptedVal =
    trspt.computeWithSecondOrderFV(velsFC,
                                   valToTransport,
                                   types,
                                   times);
  // Add contribution to RHS and LHS
  double conv(0.), contrib(0.);
  int dim(_grid->getDimension());
  for (auto cell = _grid->internalCellBegin(); cell != _grid->internalCellEnd(); cell++)
  {
    long id = cell.getId();
    int g_i = _mapGlobal.at(id);
    double diagVal = 1.;

    for (int i=0; i<dim; i++)
    {

      // Convection term
      conv = cAB_n * trsptedVal[i][id]
             + cAB_nm1 * trsptedVal[i + dim][id];

      // Prediction with Gear-scheme
      contrib = c_n * sol.VelocityCC[id][i]
                + c_nm1 * solPrev.VelocityCC[id][i]
                + c_f * (-conv - gradPressure[id][i]/_rhocc[id])+c_f*divergence[id][i]/_rhocc[id];
      switch(i)
      {
      case 0:
        lapU->addMatrixValue(g_i, g_i, diagVal);
        lapU->addRHSValue(g_i,contrib);
        break;
      case 1:
        lapV->addMatrixValue(g_i, g_i, diagVal);
        lapV->addRHSValue(g_i,contrib);
        break;
      case 2:
        //FIXME: Add contrib for third dimension
        break;
      default:
        break;
      }
    }
  }

  // Solve on each components
  // FIXME: May be an augmented matrix would be better
  this->solvePrediction(solNext.Ux, solNext.Uy);
  //FIXME: Add solve for the third dimension


  // Put Ux^{*} and Uy^{*} in VelocityCC vector.
  // Communicate ghost after. FIXME: A reflexion on useful communication
  // has to be lead.
  for (auto cell = _grid->internalCellBegin(); cell != _grid->internalCellEnd(); cell++)
  {
    const long& id = cell.getId();

    solNext.VelocityCC[id] = {solNext.Ux[id], solNext.Uy[id], 0.};
  }

  if (Prediction::_pvComm_only_one_init_dirty!=10101)
  {
    Prediction::_pvComm = new UserDataComm<double>(*_grid, 3);
    Prediction::_pvComm_only_one_init_dirty=10101;
  }

  Prediction::_pvComm->update();
  Prediction::_pvComm->communicate(solNext.VelocityCC);

  //FIXME: Add Imposed velocity if there is fluid-structure interaction

}

void Prediction::solvePrediction(PiercedVector<double>& Ux,
                            PiercedVector<double>& Uy)
{
  lapU->solveLaplacian();
  lapV->solveLaplacian();

  PiercedVector<double> solU = lapU->getSolutionPV();
  PiercedVector<double> solV = lapV->getSolutionPV();

  for (auto cell = _grid->internalCellBegin(); cell != _grid->internalCellEnd(); cell++)
    {
      const long& id = cell.getId();

      Ux[id]=solU[id];
      Uy[id]=solV[id];
    }
}
}
