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

#ifndef __COMMONTOOLS_HPP__
#define __COMMONTOOLS_HPP__

/**
 * @file CommonTools.hpp
 *
 * @brief This file contains TimeStep namespace
 *
 * @author Antoine Gerard
 * @date 2019-08-02
 *
 */

#include "Grid.hpp"
#include "SimulationParameters.hpp"
#include "NeosSolution.hpp"

namespace neos {
/**
 * @namespace TimeStep
 * @brief Namespace to define useful time-stepping tools
 */
namespace CommonTools
{
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T> >& v);

template <typename T>
T getGlobalOp(Grid* grid, const std::vector<T>& v, MPI_Op op);

double getTimeStep(Grid* grid,
                   SimulationParameters* param,
                   const double &t,
                   const std::vector<std::vector<double> > Velocity);

void update(Solution& SolPrev, Solution& Sol);

void update(Solution& SolPrev, Solution& Sol, Solution& SolNext);
}
}
#endif
