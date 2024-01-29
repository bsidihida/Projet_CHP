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
 * @file   Prediction.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Mon Sep  9 13:38:26 2019
 *
 * @brief This file contains prediction class
 *
 * @copyright Inria
 */
#ifndef __PREDICTION_HPP__
#define __PREDICTION_HPP__

#include "Grid.hpp"
//#include "LaplacianFactory2.hpp"
#include "SimulationParameters.hpp"
#include "StencilBuilder.hpp"
#include "UserDataComm.hpp"
#include "NeosSolution.hpp"


#include <common.hpp>
#include <Gradient.hpp>
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

namespace neos {

  
class Laplacian2 {

protected:
  Grid *_grid;
  std::unordered_map<long, long> _mapGlobal;
  double _positionCoef;
  Gradient* _gradient;
  StencilBuilder* _stencils;
  bool _hasStencil, _neumann, _buildMatrices;
  std::vector<double> _sol_v;
  PiercedVector<double> _lambda;
  double (*_lambdaf)(double,double,double);
  void putRHS(const PiercedVector<double> &rhs, RHSOpType);

  static UserDataComm<double>* _userComm;
  static int _userComm_only_one_init_dirty;

public:
Laplacian2() {
}
/**
 * @brief LaplacianPetsc2 class constructor
 *
 * @param[in] grid Pointer on the grid
 */
Laplacian2(Grid *grid,
          bool buildMatrices=true);

Laplacian2(Grid *grid,
          StencilBuilder* stencils);

virtual ~Laplacian2() {
}
virtual void buildMatrices(int argc, char **argv) = 0;

virtual void setMatrixValue(int row, int col, double val) = 0;
virtual void addMatrixValue(int row, int col, double val) = 0;
virtual void zeroMatrix() = 0;
virtual void setRHSValue(int row, double val) = 0;
virtual void addRHSValue(int row, double val) = 0;
virtual double getRHSValue(const int row) = 0;
virtual void zeroRHS() = 0;
virtual void syncRHS() = 0;
virtual double solveLaplacian() = 0;
virtual double solveAgainLaplacian() = 0;
virtual std::vector<double> getSolution() = 0;
virtual std::vector<double> error() = 0;
virtual void dumpMatrix() = 0;
virtual void dumpRHS() = 0;

void addSolutionToVTK(std::string tag);
void buildMatrix();
void buildFVMatrix(PiercedVector<double>& kappaCC,
                              PiercedVector<double>& kappaFC,PiercedVector<double>& rho1,
                              double t,
                              Var type);
void toggleNeumann(bool neumann);
void computeRHSHeat2D();
void setStencils(StencilBuilder* stencils);
void setLambda(std::vector<double>* lambda);
void setLambda(double (*lambdaf)(double,double,double))  ;
void setDirichletCondition(double pos=BCPosition::OnInterface);
void setDirichletConditionRHS(double t=0.0);

void setRHS(const PiercedVector<double> &rhs);
void setRHS(const std::vector<double> &rhs);
void setRHS(double (*callback)(int, double));
void addRHS(const PiercedVector<double> &rhs);
void addRHS(const std::vector<double> &rhs);
void penalize(const std::vector<double> &phi);
void penalizeAtOrder1(const PiercedVector<double> &fctInd,
                                 double coeffPen);

PiercedVector<double> getSolutionPV();
std::vector<std::array<double,3> > getGradientFromSolution();

void setRHS()
{
  this->computeRHSHeat2D();
}

/**
 * @brief Solve the Laplacian2
 */
PiercedVector<double> solve()
{
  this->buildMatrix();
  this->setDirichletCondition();
  this->solveLaplacian();
  return this->getSolutionPV();
}

/**
 * @brief For a given vlue in Cell Center, compute this value on Interface
 *
 * @brief[in] Cell Interface
 *
 * @return Interpolate value on Interface
 */
double computeLambdaOnInterface(bitpit::Interface &inter);

/**
 * @brief Set the BCPosition
 *
 * @brief[in] bcPos Enum to set the BCPosition on the interface or outer the cell
 */
// TODO: obsolete, need to remove
void setBCPositionCoef(double bcPos){
    _positionCoef = bcPos;
  }
};


class LaplacianPetsc2: public Laplacian2
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
    LaplacianPetsc2(Grid *grid,
              bool buildMatrices_=true,
              int argc=0,
              char **argv=NULL) :
              Laplacian2(grid, buildMatrices_)
              {
                if (_buildMatrices) {
                  buildMatrices(argc, argv);
                }
                //PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
                PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
              }

    LaplacianPetsc2(Grid *grid,
              StencilBuilder* stencils,
              int argc=0,
              char **argv=NULL) :
              Laplacian2(grid, stencils)
              {
                if (_buildMatrices) {
                  buildMatrices(argc, argv);
                }
                PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
              }
/**
 * @brief LaplacianPetsc2 class destructor
 */
~LaplacianPetsc2();

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
 * @brief Solve the LaplacianPetsc2
 *
 * @brief[out] return residual norm
 */
double solveLaplacian();

/**
 * @brief Solve again the LaplacianPetsc2 (same system, different RHS)
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











class LaplacianFactory2
{
public:
/**
 * @brief LaplacianFactory2 class constructor
 *
 */
LaplacianFactory2() {
}

/**
 * @brief LaplacianFactory2 class destructor
 *
 */
~LaplacianFactory2() {
}

/**
 * @brief Return the good type of Laplacian2 class
 *
 * @param[in] ltype lap_type enum to select the kind of Laplacian2 to create
 * @param[in] grid Pointer on the grid
 * @param[out] return a Laplacian2 class with the good Laplacian2 instance
 */
static Laplacian2 *get(lapType ltype, solverType lap,Grid *grid,
                       int argc, char** argv=NULL);

/**
 * @brief Return the good type of Laplacian2 class
 *
 * @param[in] ltype lap_type enum to select the kind of Laplacian2 to create
 * @param[in] grid Pointer on the grid
 * @param[out] return a Laplacian2 class with the good Laplacian2 instance
 */
//static Laplacian2 *get(lapType ltype, solverType lap, Grid *grid, int argc, char** argv);

/**
 * @brief Return the good type of Laplacian2 class
 *
 * @param[in] ltype lap_type enum to select the kind of Laplacian2 to create
 * @param[in] grid Pointer on the grid
 * @param[out] return a Laplacian2 class with the good Laplacian2 instance
 */
static Laplacian2 *get(lapType ltype, solverType lap, Grid *grid);
};






class Prediction
{
public:
/**
 * @brief Constructor
 *
 * @param[in] grid: Pointer to the current grid
 * @param[in] stencils: Pointer to the StencilBuilder used to compute
 * stencils
 * @param[in] param: Pointer to the SimulationParameters.
 */
Prediction(Grid* grid,
           StencilBuilder* stencils = NULL,
           SimulationParameters* param = NULL);

/**
 * @brief destructor
 *
 */
~Prediction();

/**
 * @brief compute the prediction step in a Navier-Stokes problem
 *        See Antoine Fondaneche's Navier Stokes report.
 * @param[in]  SolPrev: Reference to the solution structure at time n-1
 * @param[in]  Sol: Reference to the solution structure at time n
 * @param[out] SolNext: Reference to the solution structure at time n+1.
 * 
 */


void ComputePredictionStep(Solution& SolPrev,
                           Solution& Sol,
                           Solution& SolNext,PiercedVector<double>& rho1,PiercedVector<double>& mu1);

/**
 * @brief Solve the prediction linear system. LHS*res = RHS.
 *
 * @param[out] res: Solution
 * @param[in]  LHS: left hand side matrix
 * @param[in]  RHS: right hand side vector
 */
void solvePrediction(PiercedVector<double> &Ux,
                     PiercedVector<double> &Uy);


private:
Grid* _grid;                   /**< Pointer to the current grid */

std::unordered_map<long, long> _mapGlobal;     /**< local to global mapping  */

SimulationParameters* _param;     /**< Pointer to simulation parameters */
bool _hasParam;                 /**<  Does the user give parameters */
Laplacian2 *lapU;
Laplacian2 *lapV;

StencilBuilder* _stencils;      /**< Pointer to the StencilBuilder */
static UserDataComm<double>* _pvComm;     /**< MPI communicator */
static int _pvComm_only_one_init_dirty;     /**< Integer to know if the
                                               communicator was initialized */
};
}

#endif
