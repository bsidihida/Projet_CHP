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

#ifndef __SPHERE_HPP__
#define __SPHERE_HPP__

/**
 * @file Sphere.hpp
 * @brief This file contains the Circle class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <iostream>
#include <vector>

#include "Geo3D.hpp"

/**
 * @class Sphere
 * @brief Sphere geometry class
 */

namespace neos {

class Sphere : public Geo3D {
public:
/**
 * @brief Sphere class constructor
 *
 * @param[in] slice Number of slices
 * @param[in] precision Number of boundary points by slice
 * @param[in] diameter Diameter size
 * @param[in] xi The X initial position
 * @param[in] yi The Y initial position
 * @param[in] zi The Z initial position
 */
Sphere(double xi, double yi, double zi, double radius, int slice, int precision)
  : Geo3D(xi, yi, zi, slice, precision, radius) {
}

/**
 * @brief Compute boundary points
 *
 * Compute _precision boundary points.
 */
void compute();

};
}
#endif /* __SPHERE_HPP__ */
