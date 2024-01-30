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

#ifndef __GEO3D_HPP__
#define __GEO3D_HPP__

/**
 * @file Geo3D.hpp
 * @brief This file contains 3D Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-20
 * @copyright Inria
 */
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "CPGeometry.hpp"
#include "Grid.hpp"

/**
 * @class Geometry
 * @brief Geometry abstract class
 *
 * This class cannot be instantiated!
 * It's an interface for all others 3D geometry class (Sphere...)
 */
namespace neos {

class Geo3D : public CPGeometry {
public:
/**
 * @brief Geo3D class constructor
 *
 * @param[in] slice Number of slices
 * @param[in] precision Number of boundary points by slice
 * @param[in] diameter Diameter size
 * @param[in] xi The X initial position
 * @param[in] yi The Y initial position
 * @param[in] zi The Z initial position
 */
Geo3D(double xi, double yi, double zi, double radius, int slice, int precision)
  : CPGeometry(xi, yi, zi, radius, slice, precision) {
}

/**
 * @brief Geo3D class destructor
 */
virtual ~Geo3D() {
}

std::array<double, 3> getPoint(const size_t i) {
  if (i > _bdr.size()) {
    std::cout << "Index out of range" << std::endl;
    return {};
  }
  return _bdr[i];
}

std::array<double, 3> getPoint(const size_t s, const size_t p) {
  if (s * _precision + p > _bdr.size()) {
    std::cout << "Index out of range" << std::endl;
    return {};
  }
  return _bdr[s * _precision + p];
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

#endif /* __GEO3D_HPP__ */
