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
 * @file LaplacianPetsc.hpp
 * @brief This file contains florian LaplacianPetsc algo
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-12-05
 * @copyright Inria
 */

#ifndef __LAPLACIANPETSC_HPP__
#define __LAPLACIANPETSC_HPP__

#include <common.hpp>
#include <Laplacian.hpp>
#include <StencilBuilder.hpp>

#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

namespace neos {

/**
 * @class LaplacianPetsc
 * @brief LaplacianPetsc algo class
 *
 * This class contains Florian algo for LaplacianPetsc
 */
class LaplacianPetsc : public Laplacian
{
private:
KSP _ksp = NULL;
Mat _LaplacianMatrix;
Vec _RHS;
Vec _dum;
Vec _exact;
MatNullSpace _nullsp;
PetscMPIInt size;
PetscMPIInt rank;

void buildMatrices(int argc = 0, char **argv = NULL);

/**
 * @brief
 */
void exactHeat(Vec& vec);

public:
    LaplacianPetsc(Grid *grid,
              bool buildMatrices_=true,
              int argc=0,
              char **argv=NULL) :
              Laplacian(grid, buildMatrices_)
              {
                if (_buildMatrices) {
                  buildMatrices(argc, argv);
                }
                //PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
                PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
              }

    LaplacianPetsc(Grid *grid,
              StencilBuilder* stencils,
              int argc=0,
              char **argv=NULL) :
              Laplacian(grid, stencils)
              {
                if (_buildMatrices) {
                  buildMatrices(argc, argv);
                }
                PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
              }
/**
 * @brief LaplacianPetsc class destructor
 */
~LaplacianPetsc();

inline void setMatrixValue(int row, int col, double val)
{
  MatSetValue(_LaplacianMatrix, row, col, val, INSERT_VALUES);
}
inline void addMatrixValue(int row, int col, double val)
{
  MatSetValue(_LaplacianMatrix, row, col, val, ADD_VALUES);
}
inline void setRHSValue(int row, double val)
{
  VecSetValue(_RHS, row, val, INSERT_VALUES);
}
inline void addRHSValue(int row, double val)
{
  VecSetValue(_RHS, row, val, ADD_VALUES);
}
inline double getRHSValue(const int row)
{
  double temp;
  VecGetValues(_RHS, 1, &row, &temp);
  return temp;
}
inline void syncRHS()
{
  VecAssemblyBegin(_RHS);
  VecAssemblyEnd(_RHS);
}

void zeroMatrix();
void zeroRHS();

void dumpMatrix();
void dumpRHS();

/**
 * @brief Solve the LaplacianPetsc
 *
 * @brief[out] return residual norm
 */
double solveLaplacian();

/**
 * @brief Solve again the LaplacianPetsc (same system, different RHS)
 *
 * @brief[out] return residual norm
 */
double solveAgainLaplacian();

/**
 * @brief get Solution converted in std::vector
 */
std::vector<double> getSolution();

/*
 * @brief Compute the error (for test purpose)
 */
std::vector<double> error();

};
}
#endif /* __LaplacianPetscPETSC_HPP__ */
