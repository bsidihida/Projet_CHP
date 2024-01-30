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

#ifndef __AGEOMETRY_HPP__
#define __AGEOMETRY_HPP__

#include <array>
#include <vector>
#include <NGeometry.hpp>
#include <Grid.hpp>

namespace neos {

/**
 * @file AGeometry.hpp
 * @brief This file contains Analytic Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */
class AGeometry : public NGeometry {
public:
/**
 * @brief AGeometry class constructor
 *
 * @param[in] x The X initial position
 * @param[in] y The Y initial position
 * @param[in] z The Z initial position
 * @param[in] radius Radius size
 * @param[in] dim Dimension of the geometry
 */
AGeometry(double x, double y, double z, double radius, int dim)
  : NGeometry(x, y, z, radius)
{
  _dim = dim;
}

/**
 * @brief AGeometry class destructor
 */
virtual ~AGeometry() {
}

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 */
void updateCenter(double step, const std::vector<double> &u) {
  for (int i = 0; i < _dim; ++i) {
    _center[i] = _center[i] + step * u[i];
  }
}

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] id id of the cell (with bitpit)
 * @return phivalue is returned
 */
double getLevelSet(const long id) {
  return _phi[id];
}

/**
 * @brief Compute the levelset
 *
 * @param[in] grid Grid to compute the levelset
 */
virtual void computeLevelSet(Grid *grid) = 0;

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 *
 * @return Levelset of the point
 */
virtual double getLevelSet(const std::array<double, 3> &point) = 0;

/**
 * @brief Return the levelset
 *
 * @return Vector of levelset
 */
//    virtual std::vector<double> getLevelSet() = 0;

/**
 * @brief Copy the levelset
 *
 * @param[in] phistar Vector of levelset
 * @param[in] grid Current grid
 */
void setLevelSet(PiercedVector<double> &phistar, Grid *grid) {
  if (phistar.size() != _phi.size()) {
    std::cout << "PiercedVector do not have the same size" << std::endl;
  }
  /// @TODO remplacer par une copie ?
  for (auto & cell : grid->getCells()) {
    const long &id = cell.getId();
    _phi[id] = phistar[id];
  }
}
};
}

#endif /* __AGEOMETRY_HPP__ */
