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

#ifndef __LAPLACIAN_HPP__
#define __LAPLACIAN_HPP__
/**
 * @file Laplacian.hpp
 *
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 * @copyright Inria
 *
 * @author Laurent Facq
 * @date 2021-03-27
 * @copyright CNRS
 */

#include <common.hpp>
#include <UserDataComm.hpp>
#include <Gradient.hpp>
#include <StencilBuilder.hpp>
//#include <bitpit.hpp>

namespace neos {

class Laplacian {

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
Laplacian() {
}
/**
 * @brief LaplacianPetsc class constructor
 *
 * @param[in] grid Pointer on the grid
 */
Laplacian(Grid *grid,
          bool buildMatrices=true);

Laplacian(Grid *grid,
          StencilBuilder* stencils);

virtual ~Laplacian() {
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
                            PiercedVector<double>& kappaFC,
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
 * @brief Solve the Laplacian
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

}

#endif /* __LAPLACIAN_HPP__ */
