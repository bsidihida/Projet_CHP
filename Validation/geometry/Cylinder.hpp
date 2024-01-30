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

#ifndef __CYLINDER_HPP__
#define __CYLINDER_HPP__

/**
 * @file Cylinder.hpp
 * @brief This file contains the Cylinder class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <iostream>
#include <vector>

#include "Geo3D.hpp"

namespace neos {

/**
 * @class Cylinder
 * @brief Cylinder geometry class
 */
class Cylinder : public Geo3D {
private:
double _lenght;
public:
/**
 * @brief Cylinder class constructor
 *
 * @param[in] precision Number of boundary points
 * @param[in] xi The X initial position
 * @param[in] yi The Y initial position
 * @param[in] zi The Z initial position
 * @param[in] radius Radius size
 */
Cylinder(double xi, double yi, double zi, double radius, int slice, int precision, double lenght)
  : Geo3D(xi, yi, zi, radius, slice, precision) {
  _lenght = lenght;
}

/**
 * @brief Compute boundary points
 *
 * Compute _precision boundary points.
 */
void compute();
};
}
#endif /* __CYLINDER_HPP__ */
