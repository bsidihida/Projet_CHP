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

#ifndef __UTILS_HPP__
#define __UTILS_HPP__

/**
 * @file Utils.hpp
 * @brief Some conversion tools
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 * @copyright Inria
 */

#include <common.hpp>
#include <NeosPiercedVector.hpp>
#include <Grid.hpp>

namespace neos {

/**
 * @brief function to convert a PiercedVector to a std::vector
 *
 * @param[in] PV PiercedVector to convert
 * @Param[in] Grid Pointer to the Neos::Grid associated to the piercedvector
 */
template<typename T1>
std::vector<T1> PVtoV(const PiercedVector<T1> &PV, Grid *grid) {
  typename std::vector<T1> V(grid->nbCells());
  typename std::vector<T1>::iterator iV = V.begin();

  for (auto &cell : grid->getCells()) {
    if (cell.isInterior()) {
      const long &id = cell.getId();
      *iV = PV[id];
      ++iV;
    }
  }
  return V;
}

/**
 * @brief function to convert a std::vector to a PiercedVector
 *
 * @param[in] V std::vector to convert
 * @Param[in] Grid Pointer to the Neos::Grid associated to the future piercedvector
 */
template<typename T1>
PiercedVector<T1> VtoPV(const std::vector<T1> &V, Grid *grid) {
  typename neos::PiercedVector<T1> PV;
  typename std::vector<T1>::const_iterator iV = V.begin();

  PV.reserve(grid->nbCells());
  for (auto &cell : grid->getCells()) {
    if (cell.isInterior()) {
      const long &i = cell.getId();
      PV.emplace(i);
    }
  }

  for (auto &cell : grid->getCells()) {
    if (cell.isInterior()) {
      const long &id = cell.getId();
      PV[id] = *iV;
      ++iV;
    }
  }
  return PV;
}


}

#endif /* __UTILS_HPP__ */
