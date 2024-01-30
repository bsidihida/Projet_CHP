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

#ifndef _ILEVELSET_HPP_
#define _ILEVELSET_HPP_

/**
 * @file ILevelSet.hpp
 * @brief This file contains the LevelSet classes interface
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <Geometry.hpp>
#include <Grid.hpp>
#include <common.hpp>

namespace neos {

/**
 * @class ILevelSet
 * @brief Class to apply the Levelset method to the Geometry objects
 * @details This Class is for internal use for neos::LevelSet
 */
class ILevelSet {
private:
Grid                  *_grid;   /**< Reference to the Mesh class */
std::vector<Geometry*>        _geo;     /**< Vector of Geometry */
std::vector<double>           _lvst;    /**< levelset */

public:
/**
 * @brief LevelSet class constructor
 * The constructor cannot be call
 *
 * @param[in] grid reference to the Mesh class
 */
explicit ILevelSet(Grid *grid) : _grid(grid)
{
}


/**
 * @brief LevelSet class destructor
 */
~ILevelSet() {
}

/**
 * @brief Add a new Geometry in the vector
 *
 * @param[in] geo Pointer to the Geometry
 */
void addGeometry(Geometry *geo) {
  _geo.push_back(geo);
}

/**
 * @brief Return the first Geometry in the vector
 */
Geometry *getGeometry() {
  return _geo.front();
}

/**
 * @brief Return the levelset of the grid without update
 */
std::vector<double> get() {
  return _lvst;
}

/**
 * @brief Return the number of geometry
 */
int getGeoNbr() {
  return _geo.size();
}

/**
 * @brief Return the phi value of the point
 *
 * @param[in] cell Position of the point
 */
double getLevelSet(const NPoint & cell);

/**
 * @brief Update and return the levelset vector
 */
std::vector<double> getLevelSet();

void transport(Grid *grid, const PiercedVector<NPoint> &u, const double dt);

double getLevelSet(const long &id);
};
}
#endif /* __ILEVELSET_HPP__ */
