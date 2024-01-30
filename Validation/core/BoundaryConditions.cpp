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
 * @file   BoundaryConditions.cpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Tue Aug 13 13:59:34 2019
 *
 * @brief  This file contains BoundaryConditions class
 *
 *
 */

#include "BoundaryConditions.hpp"

namespace neos {

BoundaryConditions::BoundaryConditions()
{
}

void BoundaryConditions::addCondition(int border,
                                      std::string cond,
                                      std::function<double(const NPoint,
                                                           double)> f,
                                      Var type)
{
  std::pair<Var, int> varFace = std::make_pair(type, border);
  this->_BCConditions[varFace] = std::make_pair(cond, f);
}

std::map<std::pair<Var, int>,
         std::pair<std::string,
                   std::function<double(const NPoint,
                                        double)> > >& BoundaryConditions::getBCConditions()
{
  return this->_BCConditions;
}

std::string BoundaryConditions::getBCType(const int boundId,
                                          Var type)
{
  std::pair<Var, int> varFace = std::make_pair(type, boundId);
  return this->_BCConditions[varFace].first;
}

std::function<double(const NPoint,
                     double)> BoundaryConditions::getBCFunction(const int boundId,
                                                                Var type)
{
  std::pair<Var, int> varFace = std::make_pair(type, boundId);
  return this->_BCConditions[varFace].second;
}

}
