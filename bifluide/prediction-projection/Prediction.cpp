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
#include "Transport.hpp"
#include <chrono>
#include <numeric>
#include "Gradient.hpp"
//#include "LaplacianPetsc2.hpp"



//#include "Laplacian2.hpp"
#include "Utils.hpp"
#include <VTK.hpp>
#include "InterpolatorFactory.hpp"



namespace neos {

UserDataComm<double>* Laplacian2::_userComm;
int Laplacian2::_userComm_only_one_init_dirty=0;

Laplacian2::Laplacian2(Grid *grid_, bool buildMatrices_)
  : _grid(grid_),
  _mapGlobal(_grid->getCellGlobalMap()),
  _positionCoef(1),
  _gradient(NULL),
  _hasStencil(false),
  _neumann(true),
  _buildMatrices(buildMatrices_),
  _lambda({}),
  _lambdaf(NULL)
{}

Laplacian2::Laplacian2(Grid* grid, StencilBuilder* stencils)
  : _grid(grid),
  _mapGlobal(_grid->getCellGlobalMap()),
  _positionCoef(1),
  _gradient(NULL),
  _stencils(stencils),
  _hasStencil(stencils != NULL),
  _neumann(true),
  _buildMatrices(true),
  _lambda({}),
  _lambdaf(NULL)
{}

void Laplacian2::setStencils(StencilBuilder* stencils)
{
  _stencils = stencils;
  _hasStencil = stencils != NULL;
}

void Laplacian2::setLambda(std::vector<double>* lambda)
{
  _lambda = VtoPV(*lambda, _grid);

  if (Laplacian2::_userComm_only_one_init_dirty!=10101)
  {
    Laplacian2::_userComm = new UserDataComm<double>(*_grid, 1);
    Laplacian2::_userComm_only_one_init_dirty=10101;
  }

  Laplacian2::_userComm->update(_lambda.size());
  Laplacian2::_userComm->communicate(_lambda);

  _lambdaf = NULL;
}

void Laplacian2::setLambda(double (*lambdaf)(double,double,double))
{
  _lambdaf = lambdaf;
  if (_lambda.size()) _lambda.clear();
}

void Laplacian2::toggleNeumann(bool neumann) {
  _neumann = neumann;
}

void Laplacian2::putRHS(const PiercedVector<double> &rhs, RHSOpType operation)
{
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior())
    {
      if (operation == RHSOpType::ADD) {
        addRHSValue(gl, rhs[id]);
      } else {
        setRHSValue(gl, rhs[id]);
      }
    }
  }
}

void Laplacian2::setRHS(const PiercedVector<double> &rhs)
{
  putRHS(rhs, RHSOpType::SET);
}

void Laplacian2::setRHS(const std::vector<double> &rhs)
{
  PiercedVector<double> val = VtoPV(rhs, _grid);

  val.checkIntegrity();

  putRHS(val, RHSOpType::SET);
}

void Laplacian2::addRHS(const PiercedVector<double> &rhs)
{
  putRHS(rhs, RHSOpType::ADD);
}

void Laplacian2::addRHS(const std::vector<double> &rhs)
{
  PiercedVector<double> val = VtoPV(rhs, _grid);

  val.checkIntegrity();

  putRHS(val, RHSOpType::ADD);
}

void Laplacian2::setRHS(double (*callback)(int, double))
{
  PiercedVector<double> sol_pv = getSolutionPV();
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior())
    {
      double sol[1];
      sol[0] = callback(id, sol_pv[gl]);
      setRHSValue(gl, sol[0]);
    }
  }
}

void Laplacian2::penalize(const std::vector<double> &phi)
{
  PiercedVector<double> val = VtoPV(phi, _grid);

  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior() && val[id]!=0.0)
    {
      double t_val = val[id];
      addMatrixValue(gl, gl, t_val);
    }
  }
}

void Laplacian2::penalizeAtOrder1(const PiercedVector<double> &fctInd,
                                 double coeffPen)
{
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior() && fctInd[id] < 0.0)
    {
      addMatrixValue(gl, gl, coeffPen);
    }
  }
}

void Laplacian2::computeRHSHeat2D()
{
  std::array<double, 3> ic;
  double vp;

  for (auto &inter : _grid->getInterfaces())
  {
    const long &id=inter.getId();
    if (inter.isBorder())
    {
      int gl=_mapGlobal.at(inter.getOwner());
      ic=_grid->evalInterfaceCentroid(id);
      vp=-sinh(M_PI*(1-ic[0]))*sin(M_PI*ic[1])/sinh(M_PI)*_grid->evalInterfaceArea(id)/
          (_grid->evalCellSize(inter.getOwner())/2.0)/_grid->evalCellVolume(inter.getOwner());
      addRHSValue(gl, vp);
    }
  }
}

PiercedVector<double> Laplacian2::getSolutionPV() {
  return VtoPV(getSolution(),_grid);
}

void Laplacian2::buildMatrix()
{
  double vp, vm;

  if (_gradient == NULL) {
    _gradient = new Gradient(_grid->getDimension(), _grid);
  }

  for (auto &inter : _grid->getInterfaces())
  {
    const long &interfaceId = inter.getId();
    if (!inter.isBorder())
    {
      std::array<long, 2> owners = inter.getOwnerNeigh();

      int gr = _mapGlobal.at(owners[0]);
      int gl = _mapGlobal.at(owners[1]);

      std::vector<std::pair<long,double> > weights;

      // old version before Florian      weights = stencil.buildInterface(interfaceId);

      weights = _gradient->buildInterfaceNormalGradientStencil(interfaceId);

      double lambda_int = computeLambdaOnInterface(inter);

      for (std::size_t i = 0; i<weights.size(); ++i)
      {
        const long id = weights[i].first;
        int neigh = _mapGlobal.at(id);

        vp =  _grid->evalInterfaceArea(interfaceId)*weights[i].second/_grid->evalCellVolume(owners[0])*lambda_int;
        vm = -_grid->evalInterfaceArea(interfaceId)*weights[i].second/_grid->evalCellVolume(owners[1])*lambda_int;

        if (_grid->getCell(owners[0]).isInterior())
        {
          addMatrixValue(gr, neigh, vp);
        }
        if (_grid->getCell(owners[1]).isInterior())
        {
          addMatrixValue(gl, neigh, vm);
        }
      }
    } else {
      int gl = _mapGlobal.at(inter.getOwner());
      double vp = 0;
      if (_grid->getCell(inter.getOwner()).isInterior())
      {
        addMatrixValue(gl, gl, vp);
      }
    }
  }
}

void Laplacian2::buildFVMatrix(PiercedVector<double>& kappaCC,
                              PiercedVector<double>& kappaFC,PiercedVector<double>& rho1,
                              double t,
                              Var type)
{
  Gradient grad(_grid->getDimension(), _grid, _stencils);
  BoundaryConditions* BC(_grid->getBoundaryConditions());

  auto& stencils = (_stencils->getInterfaceGradientStencil());
  Stencil stencil;

  for (auto &inter: _grid->getInterfaces())
  {
    const long &interId = inter.getId();
    const long& ownerId = inter.getOwner();
    int g0 = _mapGlobal.at(ownerId);

    if (_hasStencil)
    {
      stencil = stencils[interId];
    }
    else
    {
      stencil = grad.buildFaceNormalGradientStencil(interId);
    }

    const double& area = _grid->evalInterfaceArea(interId);

    double volume = _grid->evalCellVolume(ownerId)*rho1[ownerId];

    if (inter.isBorder())
    {
      if (_grid->getCell(ownerId).isInterior())
      {
        if ( BC->getBCType(inter.getOwnerFace(), type)
        == "Dirichlet"  )        //Ghost boundary cell
        {
          for (size_t i = 0; i<stencil.neighs.size(); i++)
          {
            const long cellId = stencil.neighs[i];

            const int isOutsideCell = stencil.isOutsideCell[i];
            double value;
            const double& weight = stencil.weights.weights[i];
            if( isOutsideCell < 0 )
            {
              int boundaryFace = std::abs(isOutsideCell + 1 );
              NPoint ghostOctCenter =  _grid->evalCellCentroid(cellId);
              ghostOctCenter[boundaryFace/2] += pow(-1, boundaryFace + 1)
              * _grid->evalCellSize(cellId);

              value = -kappaCC[cellId] * kappaFC[interId] * weight
              * BC->getBCFunction(boundaryFace, type)(
                ghostOctCenter,
                t
              ) * area/volume;
              addRHSValue(g0, value);
            }
            else
            {
              value = kappaCC[cellId] * kappaFC[interId] * weight
              * area/volume;
              addMatrixValue(g0, g0, value);
            }
          }
        }
      }
    }    //ELSE HOMOGENEOUS NEUMANN: Nothing to add
    else
    {
      const long& neighId = inter.getNeigh();

      //int g1 = _mapGlobal.at(neighId);
      for (size_t i = 0; i<stencil.neighs.size(); i++)
      {
        const long cellId = stencil.neighs[i];
        volume = _grid->evalCellVolume(ownerId)*rho1[ownerId];

        const int isOutsideCell = stencil.isOutsideCell[i];
        double value;
        const double& weight = stencil.weights.weights[i];

        if( isOutsideCell < 0 )
        {
          if(BC->getBCType(std::abs(isOutsideCell + 1 ), type)
             == "Dirichlet")            //Ghost boundary cell
          {
            int boundaryFace = std::abs(isOutsideCell + 1 );
            NPoint ghostOctCenter =  _grid->evalCellCentroid(cellId);
            ghostOctCenter[boundaryFace/2] += pow(-1, boundaryFace + 1)
                                              * _grid->evalCellSize(cellId);

            value = -kappaCC[cellId] * kappaFC[interId] * weight
                    * BC->getBCFunction(boundaryFace, type)(ghostOctCenter,
                                                            t) * area/volume;

            if (_grid->getCell(ownerId).isInterior())
            {
              addRHSValue(g0,value);
            }

            if (_grid->getCell(neighId).isInterior())
            {
              int g1 = _mapGlobal.at(neighId);
              value *= volume;

              volume = _grid->evalCellVolume(neighId)*rho1[neighId];
              value  = -value/volume;

              addRHSValue(g1,value);
            }
          }
          //ELSE HOMOGENEOUS NEUMANN
        }
        else
        {
          int gn = _mapGlobal.at(cellId);
          value = kappaCC[cellId] * kappaFC[interId] * weight
                  * area/volume;

          if (_grid->getCell(ownerId).isInterior())
          {
            addMatrixValue(g0,gn,value);
          }

          if (_grid->getCell(neighId).isInterior())
          {
            int g1 = _mapGlobal.at(neighId);
            value *= volume;

            volume = _grid->evalCellVolume(neighId)*rho1[neighId];
            value  = -value/volume;

            addMatrixValue(g1,gn,value);
          }
        }
      }
    }
  }
}

std::vector<std::array<double,3> > Laplacian2::getGradientFromSolution()
{
  Gradient stencil(_grid->getDimension(), _grid);
  PiercedVector<double> sol = getSolutionPV();

  return stencil.computeLSGradient(sol);
}

void Laplacian2::setDirichletCondition(double pos)
{
  double vp;
  

  _neumann = false;
  for (auto &inter : _grid->getInterfaces())
    {
    const long &id = inter.getId();
    if (inter.isBorder())
    {
      int cellId= inter.getOwner();
      int gl = _mapGlobal.at(cellId);

      double lambda_int = computeLambdaOnInterface(inter);

      vp = -_grid->evalInterfaceArea(id) / (_grid->evalCellSize(inter.getOwner()) * pos) / _grid->evalCellVolume(inter.getOwner())*lambda_int;
      if (_grid->getCell(inter.getOwner()).isInterior())
	    {
        addMatrixValue(gl,gl,vp);
      }
    }
  }
}

void Laplacian2::setDirichletConditionRHS(double t)
{
  _neumann = false;
  for (auto &inter : _grid->getInterfaces())
  {
    if (inter.isBorder())
    {
      int gl = _mapGlobal.at(inter.getOwner());
      if (_grid->getCell(inter.getOwner()).isInterior())
      {
        double val, temp;
        temp = getRHSValue(gl);
        val = temp-2*t;
        setRHSValue(gl,val);
      }
    }
  }
}

void Laplacian2::addSolutionToVTK(std::string tag) {

  std::vector<double> v = PVtoV(getSolutionPV(),_grid);
  _grid->getVTK().addData(tag, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, _sol_v);
}

double Laplacian2::computeLambdaOnInterface(bitpit::Interface &inter) {
  double lambda_int = 1.0;

  if (_lambdaf) {

    NPoint xpos = _grid->evalInterfaceCentroid(inter.getId());
    lambda_int = (*_lambdaf)(xpos[0], xpos[1], xpos[2]);

  } else if (_lambda.size()){

     if (_gradient == NULL) {
       _gradient = new Gradient(_grid->getDimension(), _grid);
     }


    IInterpolator *interpo = InterpolatorFactory::get(interpoType::DISTANCEWEIGHTED);

    if (!inter.isBorder())
    {
      const long &interfaceId = inter.getId();

      std::array<long, 2> owners = inter.getOwnerNeigh();

      std::vector<std::pair<long,double> > weights;
      weights = _gradient->buildInterfaceNormalGradientStencil(interfaceId);

      NPoint xpos = _grid->evalInterfaceCentroid(interfaceId);

      lambda_int = (_lambda[owners[0]]+_lambda[owners[1]])/2.0;
      std::vector<NPoint> xref;
      std::vector<double> val_ref;

      for (std::size_t i = 0; i<weights.size(); ++i)
      {
        const long cellId = weights[i].first;

        xref.push_back(_grid->evalCellCentroid(cellId));
        val_ref.push_back(_lambda[cellId]);
      }

      if (weights.size()>0) {
        lambda_int = interpo->computeInterpolation(xpos, xref, val_ref);
      }
    } else {

      const long &interfaceId = inter.getId();
      const long &cellId = inter.getOwner();
      const bitpit::Cell cell = _grid->getCell(cellId);
      int norm = -1;
      double oppositeLambda = 0.0;
      int nbfaces = 1;

      NPoint interCentroid = _grid->evalInterfaceCentroid(interfaceId);
      NPoint cellCentroid = _grid->evalCellCentroid(cellId);

      NPoint oppositeCentroid;
      oppositeCentroid[NPX] = (cellCentroid[NPX]-interCentroid[NPX])+cellCentroid[NPX];
      oppositeCentroid[NPY] = (cellCentroid[NPY]-interCentroid[NPY])+cellCentroid[NPY];
      oppositeCentroid[NPZ] = (cellCentroid[NPZ]-interCentroid[NPZ])+cellCentroid[NPZ];

      if ((cellCentroid[NPY] == interCentroid[NPY]) && (cellCentroid[NPZ] == interCentroid[NPZ])) {
        norm = NPX;
      } else if ((cellCentroid[NPX] == interCentroid[NPX]) && (cellCentroid[NPZ] == interCentroid[NPZ])) {
        norm = NPY;
      } else {
        norm = NPZ;
      }

      for (auto &fid : _grid->findCellFaceNeighs(cellId)) {
        NPoint centroid = _grid->evalCellCentroid(fid);
        if ((
          (centroid[norm] < interCentroid[norm]) &&
          (centroid[norm] < oppositeCentroid[norm])
        ) || (
          (centroid[norm] > interCentroid[norm]) &&
          (centroid[norm] > oppositeCentroid[norm])
        ))
        oppositeLambda += _lambda[fid];
        nbfaces++;
      }

      oppositeLambda = oppositeLambda / nbfaces;

      if (nbfaces == 1) {
        lambda_int = _lambda[cellId] + (_lambda[cellId] - oppositeLambda)/2.0;
      } else {
        lambda_int =  _lambda[cellId] + (_lambda[cellId] - oppositeLambda)/1.5;
      }
    }
  }
  return lambda_int;
}

LaplacianPetsc2::~LaplacianPetsc2()
{
  if (_buildMatrices)
  {
    MatDestroy(&_LaplacianMatrix);
    KSPDestroy(&_ksp);
    VecDestroy(&_RHS);
    VecDestroy(&_dum);
    VecDestroy(&_exact);
  }
}

void LaplacianPetsc2::zeroRHS()
{
  VecSet(_RHS, 0.);
}
void LaplacianPetsc2::zeroMatrix()
{
  MatZeroEntries(_LaplacianMatrix);
}

double LaplacianPetsc2::solveLaplacian()
{
  PetscBool assembled;

  if (!MatAssembled(_LaplacianMatrix,&assembled)) {
    MatAssemblyBegin(_LaplacianMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_LaplacianMatrix, MAT_FINAL_ASSEMBLY);
  }

  VecDuplicate(_RHS, &_dum);

  KSPSetOperators(_ksp, _LaplacianMatrix, _LaplacianMatrix);
  KSPSetFromOptions(_ksp);
  KSPSetUp(_ksp);
  if (_neumann)
  {
    PetscBool isNull;
    MatNullSpaceTest(_nullsp, _LaplacianMatrix, &isNull);
    if (isNull) {
      MatNullSpaceRemove(_nullsp, _RHS);
      MatSetNullSpace(_LaplacianMatrix, _nullsp);
    }
  }

  return LaplacianPetsc2::solveAgainLaplacian();
}

double LaplacianPetsc2::solveAgainLaplacian()
{
  PetscInt nbIteration;
  PetscReal norm;

  VecAssemblyBegin(_RHS);
  VecAssemblyEnd(_RHS);

  KSPSolve(_ksp, _RHS, _dum);
  KSPGetIterationNumber(_ksp, &nbIteration);
  KSPGetResidualNorm(_ksp, &norm);
  if (nbIteration == 100000)
  {
    std::cerr<<"ERROR: LaplacianPetsc2 solver did not converged"<<std::endl;
    exit(0);
  }
  else if (nbIteration == 0)
    {
      std::cerr<<"ERROR: LaplacianPetsc2 solver did 0 iteration" <<std::endl;
      exit(0);
    }
  return (double)norm;
}

std::vector<double> LaplacianPetsc2::getSolution()
{
  _sol_v.clear();
  _sol_v.reserve(_grid->nbCells());
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    PetscInt gl = _mapGlobal.at(id);
    if (_grid->getCell(id).isInterior())
    {
      double sol;
      VecGetValues(_dum, 1, &gl, &sol);
      _sol_v.push_back(sol);
    }
  }
  return _sol_v;
}

void LaplacianPetsc2::exactHeat(Vec& vec)
{
  std::array<double, 3> ic;
  double vp;

  for (auto cell = _grid->internalCellBegin(); cell != _grid->internalCellEnd(); ++cell)
  {
    const long &id = cell.getId();
    PetscInt gl = _mapGlobal.at(id);

    ic = _grid->evalCellCentroid(id);
    vp = sinh(M_PI * (1-ic[0])) * sin(M_PI * ic[1]) / sinh(M_PI);
    VecSetValues(vec, 1, &gl, &vp, INSERT_VALUES);
  }
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);
}

std::vector<double> LaplacianPetsc2::error()
{
  exactHeat(_exact);
  double eval, sol;
  std::vector<double> err(3,0.0);

  for (auto cell = _grid->internalCellBegin(); cell!=_grid->internalCellEnd(); ++cell)
  {
    const long &id = cell.getId();
    PetscInt gl = _mapGlobal.at(id);

    VecGetValues(_dum, 1, &gl, &eval);
    VecGetValues(_exact, 1, &gl, &sol);
    err[0] += fabs(eval-sol) * _grid->evalCellVolume(id);
    err[1] += (eval-sol) * (eval-sol) * _grid->evalCellVolume(id);
    err[2] = std::max(err[2], fabs(eval-sol));
  }
  return err;
}

void LaplacianPetsc2::buildMatrices(int argc, char **argv) {

  _buildMatrices = true;

  PetscInitialize(&argc,&argv,(char *)0,"");

  VecCreate(PETSC_COMM_WORLD, &_RHS);
  VecSetType(_RHS, VECSTANDARD);

  MatCreate(PETSC_COMM_WORLD, &_LaplacianMatrix);
  MatSetType(_LaplacianMatrix, MATMPIAIJ);

  VecSetSizes(_RHS, _grid->nbCells(), PETSC_DECIDE);
  VecDuplicate(_RHS, &_dum);
  VecDuplicate(_RHS, &_exact);

  MatSetSizes(_LaplacianMatrix, _grid->nbCells(), _grid->nbCells(),
              PETSC_DECIDE, PETSC_DECIDE);
  MatMPIAIJSetPreallocation(_LaplacianMatrix, 90, PETSC_NULL, 90, PETSC_NULL);
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &_nullsp);

  MatZeroEntries(_LaplacianMatrix);
  VecZeroEntries(_RHS);

  KSPCreate(PETSC_COMM_WORLD, &_ksp);
  KSPSetType(_ksp, KSPBCGS);
  KSPSetTolerances(_ksp, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
}

void LaplacianPetsc2::dumpMatrix() {
  PetscBool assembled;

  if (!MatAssembled(_LaplacianMatrix,&assembled)) {
    MatAssemblyBegin(_LaplacianMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_LaplacianMatrix, MAT_FINAL_ASSEMBLY);
  }
  MatView(_LaplacianMatrix,PETSC_VIEWER_STDOUT_WORLD);
}

void LaplacianPetsc2::dumpRHS() {
  VecAssemblyBegin(_RHS);
  VecAssemblyEnd(_RHS);
  VecView(_RHS,PETSC_VIEWER_STDOUT_WORLD);
}



Laplacian2 *LaplacianFactory2::get(lapType ltype, solverType solver,
                                  Grid *grid,
                                  int argc, char **argv)
{
  Laplacian2 *lap = NULL;
  switch (ltype)
  {
  case lapType::FINITEVOLUME:
    lap =  new LaplacianPetsc2(grid, true, argc,argv);
    break;
  case lapType::FINITEDIFF:
    lap =  new LaplacianPetsc2(grid, true, argc,argv);
    break;
  default:
    lap =  new LaplacianPetsc2(grid, true, argc,argv);
    break;
  }
  return lap;
}


Laplacian2 *LaplacianFactory2::get(lapType ltype, solverType solver, Grid *grid)
{
  Laplacian2 *lap = NULL;
  switch (ltype)
  {
  case lapType::FINITEVOLUME:
    lap =  new LaplacianPetsc2(grid,true, 0,NULL);
    break;
  case lapType::FINITEDIFF:
    lap =  new LaplacianPetsc2(grid,true, 0, NULL);
    break;
  default:
    lap =  new LaplacianPetsc2(grid,true, 0,NULL);
    break;
  }
  return lap;
}






























  

UserDataComm<double>* Prediction::_pvComm = NULL;
int Prediction::_pvComm_only_one_init_dirty=0;

Prediction::Prediction(Grid* grid,
                       StencilBuilder* stencils,
                       SimulationParameters* param ) : _grid(grid),

  _param(param),
  _hasParam(param != NULL),
  _stencils(stencils)
{
  lapU = LaplacianFactory2::get(lapType::FINITEVOLUME, solverType::PETSC, grid);
  lapV = LaplacianFactory2::get(lapType::FINITEVOLUME, solverType::PETSC, grid);
  _mapGlobal = _grid->getCellGlobalMap();

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
                                       Solution& solNext,PiercedVector<double>& rho1,PiercedVector<double>& mu1)
{
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
  // Here we have a Laplacian2 for X and Y composant.
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
    kappaCC[cell.getId()] = 1;//mu1[cell.getId()]/rho1[cell.getId()];//
  }
  for (auto& inter: _grid->getInterfaces())
  {
    kappaFC.emplace(inter.getId());
    kappaFC[inter.getId()] = mu1[inter.getId()];
  }

  // FIXME: Separate build of matrix and vector.
  // The matrix could be computed once at each mesh refinment instead of at
  // each iteration (if there is no penalization).
  lapU->zeroMatrix();
  lapU->zeroRHS();
  lapU->buildFVMatrix(kappaCC,
                     kappaFC,rho1,
                     solNext.t + solNext.dt,
                     Var::Ux);
  lapV->zeroMatrix();
  lapV->zeroRHS();
  lapV->buildFVMatrix(kappaCC,
                     kappaFC,rho1,
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
  double conv(0.), contrib(0.),diff(0.);
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
                + c_f * (-conv - gradPressure[id][i]/rho1[id]);
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
