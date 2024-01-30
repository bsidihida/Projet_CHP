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

#ifndef __TRANSPORT_HPP__
#define __TRANSPORT_HPP__

/**
 * @file Transport.hpp
 * @brief This file contains the Transport class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 *
 * @copyright Inria
 */

#include <Grid.hpp>
#include <IInterpolator.hpp>
#include <RBF.hpp>
#include <UserDataComm.hpp>
#include <NeosPiercedVector.hpp>
#include <StencilBuilder.hpp>

namespace neos {

/**
 * @class Transport
 * @brief Transport class
 *
 * Class used to apply a transport on a vector
 */
class Transport {
private:
Grid *_grid;
enum interpoType _interpType;
double _epsilon;
IInterpolator *_interpo;
static UserDataComm<double>* _pvComm;
static int _pvComm_only_one_init_dirty;
StencilBuilder* _stencil;

double computeInterpo(PiercedVector<double> &data, const long &id, NPoint x_old);

public:
/**
 * @brief Transport class constructor
 *
 * @param[in] grid Pointeur to the grid
 * @Param[in] interpType Interpolation type
 */
Transport(Grid *grid,  enum interpoType interpType = RBF) : _grid(grid), _interpType(interpType)
{
  _epsilon = 0.1;
  _interpo = NULL;
  _stencil = nullptr;
}

/**
 * Transport class constructor
 *
 * @param[in] grid Pointer to the grid
 * @param[in] stencil Pointer to the stencilBuilder
 */
Transport(Grid *grid,  StencilBuilder* stencil);

/**
 * @brief Change the type of interpolation
 *
 * @Paran[in] interpType Interpolation type
 */
void setInterpType(enum interpoType interpType) {
  _interpType = interpType;
}

/**
 * @brief Get the type of interpolation
 *
 * @return Interpolation type
 */
enum interpoType getInterpType() {
  return _interpType;
}

/**
 * @brief Change the epsilon value (for RBF interpolation)
 *
 * @param[in] epsilon Epsilon value
 */
void setEpsilon(double epsilon) {
  _epsilon = epsilon;
}

/**
 * @brief Get the epsilon value
 *
 * @return Epsilon value
 */
double getEpsilon() {
  return _epsilon;
}

/**
 * @brief Compute the interpolation
 *
 * @Param[in] data Vector to interpolate
 * @param[in] u Speed vector
 * @param[in] dt Step size
 */
void compute(PiercedVector<double> &data,
             const std::vector<double> &u,
             const double dt);
void compute(PiercedVector<double> &data,
             const PiercedVector<NPoint> &u,
             const double dt);

/**
 * @brief  transport several vector using an order two FV method
 *
 * @param[in] velFC: Face centered velocity
 * @param[in] u: Reference on vectors to transport
 * @param[in] t: Current time
 * @param[in] type: Type of values to transport: LevSet or Values.
 *
 * @return Transported vectors
 */
std::vector<PiercedVector<double> > computeWithSecondOrderFV(
  const std::vector<PiercedVector<double> > &velFC,
  std::vector<PiercedVector<double> > &u,
  const std::vector<Var>& types,
  const std::vector<double>& times = {});

// FIXME: There is surely a best option to implement this (template ? )

/**
 * @brief transport a unique vector using an order two FV method
 *
 * @param velFC: Face centered velocity
 * @param u: Reference on vector to transport
 * @param t: Current time
 * @param type: Type of values to transport: LevSet or Values.
 *
 * @return Transported vector
 */
PiercedVector<double> computeWithSecondOrderFV(
  const PiercedVector<double> &velFC,
  PiercedVector<double> &u,
  double t=0.,
  Var type = Var::LS);
};
}

#endif /* __TRANSPORT_HPP__ */
