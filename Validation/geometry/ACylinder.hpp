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

#ifndef __ACYLINDER_HPP__
#define __ACYLINDER_HPP__

#include <array>
#include <vector>
//#include "bitpit.hpp"
#include "Grid.hpp"
#include "AGeometry.hpp"

namespace neos {

/**
 * @file ACylinder.hpp
 * @brief This file contains Analytic Cylinder class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */
class ACylinder : public AGeometry {
private:
int _axe;   /**< Axe of the cylinder */

public:
/**
 * @brief ACylinder class constructor
 *
 * @param[in] x The X initial position
 * @param[in] y The Y initial position
 * @param[in] z The Z initial position
 * @param[in] radius Radius size
 * @param[in] dim Dimension of the geometry
 * @param[in] axe Axe of the cylinder
 */
ACylinder(double x, double y, double z, double radius, int dim, int axe)
  : AGeometry(x, y, z, radius, dim) {
  _axe = axe;
}

using AGeometry::getLevelSet;

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 *
 * @return Levelset of the point
 */
double getLevelSet(const NPoint &cell);

/**
 * @brief Compute the levelset
 *
 * @param[in] grid Grid to compute the levelset
 */
void computeLevelSet(Grid *grid);

/**
 * @brief Set the axe
 *
 * @param[in] axe Cylinder new axe
 */
void setAxe(int axe) {
  _axe = axe;
}
};
}

#endif /* __ACYLINDER_HPP__ */
