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
 * @file LaplacianFactory.hpp
 * @brief This file contains a Factory design pattern to return
 * the good laplacian type
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-11-23
 * @copyright Inria
 */

#ifndef __LAPLACIANFACTORY_HPP_
#define __LAPLACIANFACTORY_HPP_

#include <common.hpp>
#include <Laplacian.hpp>
#include <Grid.hpp>

namespace neos {

/**
 * @class LaplacianFactory
 * @brief Laplacian algo class
 *
 * This class contains Florian algo for laplacian
 */
class LaplacianFactory
{
public:
/**
 * @brief LaplacianFactory class constructor
 *
 */
LaplacianFactory() {
}

/**
 * @brief LaplacianFactory class destructor
 *
 */
~LaplacianFactory() {
}

/**
 * @brief Return the good type of laplacian class
 *
 * @param[in] ltype lap_type enum to select the kind of laplacian to create
 * @param[in] grid Pointer on the grid
 * @param[out] return a Laplacian class with the good laplacian instance
 */
static Laplacian *get(lapType ltype, solverType lap,Grid *grid,
                       int argc, char** argv=NULL);

/**
 * @brief Return the good type of laplacian class
 *
 * @param[in] ltype lap_type enum to select the kind of laplacian to create
 * @param[in] grid Pointer on the grid
 * @param[out] return a Laplacian class with the good laplacian instance
 */
//static Laplacian *get(lapType ltype, solverType lap, Grid *grid, int argc, char** argv);

/**
 * @brief Return the good type of laplacian class
 *
 * @param[in] ltype lap_type enum to select the kind of laplacian to create
 * @param[in] grid Pointer on the grid
 * @param[out] return a Laplacian class with the good laplacian instance
 */
static Laplacian *get(lapType ltype, solverType lap, Grid *grid);
};
}

#endif /*  __LAPLACIANFACTORY_HPP_ */
