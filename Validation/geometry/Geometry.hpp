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

#ifndef __GEOMETRY_HPP__
#define __GEOMETRY_HPP__

/**
 * @file Geometry.hpp
 * @brief This file contains Geometry abstract class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <vector>
#include <array>
#include <common.hpp>
#include <Grid.hpp>
#include <NeosPiercedVector.hpp>

/**
 * @class Geometry
 * @brief Geometry abstract class
 *
 * This class cannot be instantiated!
 * It's an interface for all others geometry class
 *
 * getLevelSet functions are to be implemented in child class
 */

namespace neos {

class Geometry {
public:
/**
 * @brief Geometry class constructor
 *
 * The constructor cannot be call
 */
Geometry() {
}

/**
 * @brief Geometry class destructor
 */
virtual ~Geometry() {
}

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 */
virtual double getLevelSet(const std::array<double, 3> &point) = 0;
virtual double getLevelSet(const long id) = 0;
virtual void computeLevelSet(Grid *grid) = 0;
virtual PiercedVector<double> &getPhi() = 0;
/**
 * @brief Return the levelset
 *
 * @param[in] slice Number of slices
 */
//    virtual std::vector<double> getLevelSet() = 0;

};
}
#endif /* __GEOMETRY_HPP__ */
