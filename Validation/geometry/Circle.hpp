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

#ifndef __CIRCLE_HPP__
#define __CIRCLE_HPP__

/**
 * @file Circle.hpp
 * @brief This file contains the Circle class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <iostream>
#include <vector>

#include "Geo2D.hpp"

namespace neos {

/**
 * @class Circle
 * @brief Circle geometry class
 */
class Circle : public Geo2D {
public:
/**
 * @brief Circle class constructor
 *
 * @param[in] precision Number of boundary points
 * @param[in] xi The X initial position
 * @param[in] yi The Y initial position
 * @param[in] radius Radius size
 */
Circle(double xi, double yi, double radius, int precision)
  : Geo2D(xi, yi, radius, precision) {
}

/**
 * @brief Compute boundary points
 *
 * Compute _precision boundary points.
 */
void compute();
};
}

#endif /* __CIRCLE_HPP__ */
