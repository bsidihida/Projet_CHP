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

#ifndef __GEO2D_HPP__
#define __GEO2D_HPP__

/**
 * @file Geo2D.hpp
 * @brief This file contains 2D Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-20
 * @copyright Inria
 */

#include <iostream>
#include <stdlib.h>
#include <vector>
#include "CPGeometry.hpp"
#include "Grid.hpp"

/**
 * @class Geo2D
 * @brief Geo2D abstract class
 *
 * This class cannot be instantiated!
 * It's an interface for all others 2D geometry class (Circle...)
 */
namespace neos {

class Geo2D : public CPGeometry {
public:
Geo2D(double xi, double yi, double radius, int precision)
  : CPGeometry(xi, yi, radius, precision) {
}

/**
 * @brief Geo2D class destructor
 */
virtual ~Geo2D() {
}

std::array<double, 3> getPoint(const size_t i) {
  if (i > _bdr.size()) {
    std::cout << "Index out of range" << std::endl;
    return {};
  }
  return _bdr[i];
}

/**
 * @brief Compute the levelset
 *
 * @param[in] grid Grid to compute the levelset
 */
void levelSet(Grid *grid);

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 */
double getLevelSet(const std::array<double, 3> &cell);
};

}
#endif /* __GEO2D_HPP__ */
