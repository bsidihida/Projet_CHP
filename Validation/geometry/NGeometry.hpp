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

#ifndef __NGEOMETRY_HPP__
#define __NGEOMETRY_HPP__

/**
 * @file NGeometry.hpp
 * @brief This file contains Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <iostream>
#include <vector>
#include <functional>
#include <Geometry.hpp>

/**
 * @class NGeometry
 * @brief NGeometry abstract class
 *
 * This class cannot be instantiated!
 * It's an interface for all others neos geometry class
 */

namespace neos {

class NGeometry : public Geometry {
protected:
std::array<double,3>          _initialCenter;   /**< The initial X, Y, Z position */
std::array<double,3>          _center;          /**< The current X, Y, Z position */
double _radius;                                         /**< Radius size */
int _dim;                                               /**< Geometry dimension */
PiercedVector<double> _phi;             /**< Levelset vector */

public:
/**
 * @brief Geometry class constructor
 *
 * The constructor cannot be call
 *
 * @param[in] x The X initial position
 * @param[in] y The Y initial position
 * @param[in] radius Radius size
 */
NGeometry(double x, double y, double radius) : _initialCenter({ { x, y, 0 } }), _center({ { x, y, 0 } }), _radius(radius) {
  _dim       = 2;
}

/**
 * @brief Geometry class constructor
 *
 * The constructor cannot be call
 *
 * @param[in] x The X initial position
 * @param[in] y The Y initial position
 * @param[in] z The Z initial position
 * @param[in] radius Radius size
 */
NGeometry(double x, double y, double z, double radius) : _initialCenter({ { x, y, z } }), _center({ { x, y, z } }), _radius(radius) {
  _dim       = 3;
}

/**
 * @brief Geometry class destructor
 */
virtual ~NGeometry() {
}

/**
 * @brief Return the center coordinate
 *
 */
std::array<double, 3> getCenter() {
  return _center;
}

/**
 * @brief Return the geometry dimension
 *
 */
int getDim() {
  return _dim;
}

/**
 * @brief Return the center initial coordinate
 *
 */
std::array<double, 3> getInitialCenter() {
  return _initialCenter;
}

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 */
virtual double getLevelSet(const std::array<double, 3> &point) = 0;


/**
 * @brief Return the levelset of cell
 *
 * @param[in] id id of the cell
 *
 * @return levelset value
 */
double getLevelSet(const long id) {
  return _phi[id];
}

/// @TODO renommer en GetLevelSet ?
PiercedVector<double> &getPhi() {
  return _phi;
}
};
}
#endif /* __NGEOMETRY_HPP__ */
