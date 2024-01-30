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

#ifndef __BOUNDARYCONDITIONS_HPP__
#define __BOUNDARYCONDITIONS_HPP__

/**
 * @file   BoundaryConditions.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Tue Aug 13 13:55:17 2019
 *
 * @brief  This file contains BoundaryConditions class
 *
 *
 */

#include <unordered_map>
#include <iostream>
#include <functional>
#include <array>
#include <vector>
#include <map>

#include <common.hpp>

namespace neos {

/**
 *@class BoundaryConditions
 *@brief Class to handle boundary conditions
 *
 *
 **/

class BoundaryConditions
{
public:
/**
 * Default constructor
 *
 */
BoundaryConditions();

/**
 * Destructor
 *
 */
~BoundaryConditions(){
  ;
}

/**
 * @brief Add a boundary condition on border
 *
 * @param[in] border Local number of the border
 * @param[in] cond "Dirichlet" or "Neumann"
 * @param[in] f Function used as boundary conditions
 */
// FIXME: Only work for rectangular domain with a condition on each side.
void addCondition(int border,
                  std::string cond,
                  std::function<double(const NPoint, double)> f,
                  Var type = Var::LS);

/**
 * Get boundary conditions function
 *
 *
 * @return Boundary conditions in form <Condition type, function>
 */
std::map<std::pair<Var, int>,
         std::pair<std::string,
                   std::function<double(const NPoint,
                                        double)> > >& getBCConditions();

std::string getBCType(const int boundId,
                      Var type = Var::LS);

std::function<double(const NPoint,
                     double)> getBCFunction(const int boundId,
                                            Var type = Var::LS );


private:
std::map<std::pair<Var, int>,
         std::pair<std::string,
                   std::function<double(const NPoint,
                                        double)> > > _BCConditions;
/**< map which store BC type and functions for each side of the domain*/
};
}

#endif
