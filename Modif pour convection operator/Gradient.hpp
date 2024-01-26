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

#ifndef __GRADIENT_HPP__
#define __GRADIENT_HPP__

/**
 * @file Gradient.hpp
 *
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 * @copyright Inria
 */

//#include "bitpit.hpp"
#include <common.hpp>

#include <InterpolatorFactory.hpp>
#include <UserDataComm.hpp>
#include <Grid.hpp>
#include <tuple>
#include <NeosPiercedVector.hpp>
#include <StencilBuilder.hpp>

namespace neos {

class Gradient {
private:
Grid *_grid;
StencilBuilder *_stencil;
const bool _hasStencil;

protected:
std::unordered_map<long, long> _mapGlobal;
int _dim;
interpoType _iType = interpoType::DISTANCEWEIGHTED;
static int _pvComm_only_one_init_dirty;
static UserDataComm<double>* _userComm;
static int _userComm_only_one_init_dirty;

public:
// Gradient(bitpit::VolumeKernel *grid, interpo_type iType) {
Gradient(int dim, Grid *grid, StencilBuilder* stencil = nullptr);
~Gradient() {
}
std::vector<std::pair<long,double> > buildVertexStencil(const long ownerId, const int localVertex);
std::vector<std::array<double,3> >   computeFVGradient(PiercedVector<double> &val);
std::vector<std::pair<long,double> > computeNormalGradient(const long interfaceId, std::vector<std::vector<std::pair<long,double> > > &grad);
std::vector<std::vector<std::pair<long,double> > > buildInterfaceGradientStencil(const long interfaceId);
std::vector<std::pair<long,double> > buildInterfaceNormalGradientStencil(const long interfaceId);
std::map<std::pair<int,int>,double> buildSizeMap(const long interfaceId);

/**
 * @brief compute a cell centered gradient using a least-square
 * interpolation (for a unique function)
 * @param[in] val cell-centered value of the function we want to compute
 *             gradient
 * @return Gradient as a Vector
 */
std::vector<std::array<double,3> > computeLSGradient(PiercedVector<double> &val);

/**
 * @brief compute a cell centered gradient using a least-square
 * interpolation (for a unique function)
 * @param[in] val cell-centered value of the function we want to compute
 *             gradient
 * @return Gradient as a Vector
 */
std::vector<std::array<double,3> > computeLSGradient(const std::vector<double> &val);

std::array<double,3> computeLSCellGradient(long id, PiercedVector<double> &data);
void buildLSCellGridMatrix(long id, std::vector<long> &neigh,double* A);
void buildLSCellRHS(long id, std::vector<long>
                    &neigh,PiercedVector<double> &data, double *F);

/**
 * @brief compute a cell centered gradient using a least-square
 * interpolation (for a unique function)
 *
 * @param val: cell-centered value of the function we want to compute
 *             gradient
 * @param t : time at which the gradient is computed
 * @param type: type of function for boundary condition
 *
 * @return Gradient as a pierced vector
 */
PiercedVector<std::array<double, 3> > computeCCLSGradient(
  PiercedVector<double>& val,
  const double t = 0.,
  const Var type = Var::LS);

/**
 * @brief compute a cell centered gradient using a least-square
 * interpolation (for several function)
 *
 * @param val: cell-centered values of the functions we want to compute
 *             gradient
 * @param types: Vector of types for each function
 * @param times: Vector of times corresponding to each current time for each
 *               function
 *
 * @return A vector of each gradient as PiercedVector
 */
std::vector<PiercedVector<std::array<double, 3> > > computeCCLSGradient(
  std::vector<PiercedVector<double> >& val,
  const std::vector<Var>& types,
  const std::vector<double>& times);

/**
 * @brief Build a stencil to compute the gradient at cell cellId
 *
 * @param cellId: index of the cell we want to compute stencil
 *
 * @return The Stencil struct for the cell cellId
 */
Stencil buildCellGradientStencil(const long& cellId);


/**
 * @brief Compute a face centered normal gradient using a least-square
 *        interpolation
 *
 * @param val: cell-centered value of the function we want to compute FC
 * normal gradient.
 * @param t: time at which the gradient is computed
 * @param type: type of the function for boundary conditions
 *
 * @return a PiercedVector containing the Face centered normal gradient
 */
PiercedVector<double> computeFCLSGradient(PiercedVector<double>& val,
                                          const double t,
                                          const Var type);

/**
 * @brief Compute a face centered normal gradient using a least-square
 *        interpolation (for several functions)
 *
 * @param val: vector of cell-centered values of functions we want to
 * compute FC normal gradient.
 * @param types: vector of types for boundary conditions of each function
 * @param times: vector of times at which each gradient is computed.
 *
 * @return
 */
std::vector<PiercedVector<double> > computeFCLSGradient(
  std::vector<PiercedVector<double> >& val,
  const std::vector<Var>& types,
  const std::vector<double>& times);

/**
 * @brief Build a stencil to compute the normal gradient at interface
 * interId
 *
 * @param interId: index of the interface we want to compute stencil
 *
 * @return The Stencil struct for the interface interId
 */
Stencil buildFaceNormalGradientStencil(const long& interId);

/**
 * @brief Compute a finite volume cell-centered gradient.
 *        FIXME: not implemented yet for non-uniform grid
 *
 * @param val: cell-centered value of the function we want to compute CC
 *             gradient.
 * @param t: time at which the gradient is computed
 * @param type: type for boundary conditions of the function. x
 *
 * @return
 */
PiercedVector<std::array<double, 3> > computeFVGradient(
  PiercedVector<double>&  val,
  const double t,
  const Var type = Var::LS);


static inline double DotProduct(const std::array<double,3> &u, const std::array<double,3> &v){
  //int size=std::min(u.size(),v.size());
  double dotproduct=0.0;

  for (int i=0; i<3; ++i)
    dotproduct+=u[i]*v[i];

  return dotproduct;
}



};
}

#endif /* __GRADIENT_HPP__ */
