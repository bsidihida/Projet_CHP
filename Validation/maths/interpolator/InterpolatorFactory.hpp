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

/**
 * @file InterpolatorFactory.hpp
 * @brief This file contains a Factory design pattern to return
 * the good interpolator algorithm
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-11-23
 * @copyright Inria
 */

#ifndef __INTERPOFACTORY_HPP__
#define __INTERPOFACTORY_HPP__

#include "common.hpp"
#include "IInterpolator.hpp"

namespace neos {

/**
 * @class InterpolatorFactory
 * @brief Interpolator algo class
 *
 * This class contains Florian algo for laplacian
 */
class InterpolatorFactory
{
public:
/**
 * @brief InterpolatorFactory class constructor
 *
 */
InterpolatorFactory() {
}

/**
 * @brief InterpolatorFactory class destructor
 *
 */
~InterpolatorFactory() {
}

/**
 * @brief Return the good type of laplacian class
 *
 * @param[in] itype interpoType enum to select the kind of interpolator class to create
 * @param[out] return a IInterpolator class with the good laplacian instance
 */
static IInterpolator *get(interpoType itype);
};
}

#endif /*  __INTERPOFACTORY_HPP__ */
