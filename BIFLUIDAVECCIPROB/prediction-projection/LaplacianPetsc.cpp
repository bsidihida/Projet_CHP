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
 * @file LaplacianPetsc.cpp
 * @brief This file contains florian laplacian algo
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-12-05
 * @copyright Inria
 */

#include "LaplacianPetsc.hpp"
#include "Gradient.hpp"
#include "Utils.hpp"
#include <numeric>

namespace neos {

LaplacianPetsc::~LaplacianPetsc()
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

void LaplacianPetsc::zeroRHS()
{
  VecSet(_RHS, 0.);
}
void LaplacianPetsc::zeroMatrix()
{
  MatZeroEntries(_LaplacianMatrix);
}

double LaplacianPetsc::solveLaplacian()
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

  return LaplacianPetsc::solveAgainLaplacian();
}

double LaplacianPetsc::solveAgainLaplacian()
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
    std::cerr<<"ERROR: LaplacianPetsc solver did not converged"<<std::endl;
    exit(0);
  }
  else if (nbIteration == 0)
    {
      std::cerr<<"ERROR: LaplacianPetsc solver did 0 iteration" <<std::endl;
      exit(0);
    }
  return (double)norm;
}

std::vector<double> LaplacianPetsc::getSolution()
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

void LaplacianPetsc::exactHeat(Vec& vec)
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

std::vector<double> LaplacianPetsc::error()
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

void LaplacianPetsc::buildMatrices(int argc, char **argv) {

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

void LaplacianPetsc::dumpMatrix() {
  PetscBool assembled;

  if (!MatAssembled(_LaplacianMatrix,&assembled)) {
    MatAssemblyBegin(_LaplacianMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_LaplacianMatrix, MAT_FINAL_ASSEMBLY);
  }
  MatView(_LaplacianMatrix,PETSC_VIEWER_STDOUT_WORLD);
}

void LaplacianPetsc::dumpRHS() {
  VecAssemblyBegin(_RHS);
  VecAssemblyEnd(_RHS);
  VecView(_RHS,PETSC_VIEWER_STDOUT_WORLD);
}

}
