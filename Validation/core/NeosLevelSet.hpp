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

#ifndef __LEVELSET_HPP__
#define __LEVELSET_HPP__

/**
 * @file LevelSet.hpp
 * @brief This file contains the LevelSet meta classe
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <map>
#include <atomic>
#include <vector>
#include <common.hpp>
#include <Geometry.hpp>
#include <ILevelSet.hpp>
#include <NeosPiercedVector.hpp>


namespace neos {

/**
 * Return an unique id
 */
class Utils {
public:
/**
 * @brief Return an unique id
 */
static int getUid() {
  static std::atomic<std::uint32_t> uid { 0 };
  return ++uid;
}
};


/**
 * @class LevelSet
 * @brief Class to apply the Levelset method to the Geometry objects
 * @details Class Levelset store vector of neos::ILevelSets to get the LevelSet of a selection
 * of neos::Geometry objetcs (neos::ILevelSets calculate Levelset of each neos::Geometry stored)
 * you can store neos::Geometry with a tag to name geometries, and each neos::Geometry has an unique uid
 */
class LevelSet {
private:
Grid                  *_grid;   /**< Reference to the Mesh class */
std::map<int, ILevelSet*>   _lvlst;     /**< Vector of ILevelSet */
std::map<std::string, int>    _lvlstMap;   /**< association map of TAG and id */
int nbGeo;

public:
/**
 * @brief LevelSet class constructor
 *
 * @param[in] grid reference to the Mesh class
 */
explicit LevelSet(Grid *grid) : _grid(grid)
{
  nbGeo = 0;
}

/**
 * @brief LevelSet class destructor
 */
~LevelSet()
{
  for (std::map<int, ILevelSet*>::iterator itr = _lvlst.begin(); itr != _lvlst.end(); ++itr) {
    delete itr->second;
  }
}


/**
 * @brief Return the current nb of Geometry
 */
int countGeometry() {
  return nbGeo;
}

/**
 * @brief Add a Geometry
 *
 * @param[in] geo Pointer to the Geometry
 * @return Unique Id of geometry stored in Levelset
 */
int addGeometry(Geometry *geo);

/**
 * @brief Add a Geometry
 *
 * @param[in] geo Pointer to the Geometry
 * @param[in] tag Tag of the geometry
 * @return Unique Id of geometry stored in Levelset
 */
int addGeometry(Geometry *geo, std::string tag);

/**
 * @brief Delete a levelset
 *
 * @param[in] id Unique Geometry id
 */
void delGeometry(int id);

/**
 * @brief Delete a group of levelset
 *
 * @param[in] tag Group tag
 */
void delGeometry(std::string tag);

/**
 * @brief Return a Geometry
 *
 * @param[in] idGeo Unique id of the geometry
 */
Geometry *getGeometry(int idGeo);

/**
 * @brief Return the id of a given tag
 *
 * @param[in] tag Geometry tag
 */
int getId(std::string tag);

/**
 * Print the Geometry map (for debug)
 */
void printMaps();

/**
 * @brief Return the levelset value of an given point
 *
 * @param[in] cell Coordinate of the point
 */
double getLevelSet(const NPoint &cell);

/**
 * @brief Return the levelset value of an given point for a geometry
 *
 * @param[in] cell Coordinate of the point
 * @param[in] id Unique id of the Geometry
 */
double getLevelSet(const NPoint &cell, int id);

/**
 * @brief Return the levelset value of an given point for a tag
 *
 * @param[in] cell Coordinate of the point
 * @param[in] tag Tag of the set of geometry
 */
double getLevelSet(const NPoint &cell, std::string tag);

/**
 * @brief Return the levelset of the grid
 */
std::vector<double> getLevelSet();

/**
 * @brief Return the levelset of the grid for a given Geometry
 *
 * @param[in] id Unique id of the Geometry
 */
std::vector<double> getLevelSet(int id);

/**
 * @brief Return the levelset of the grid for a given set of Geometry
 *
 * @param[in] tag Tag of the Geometrie's set
 */
std::vector<double> getLevelSet(std::string tag);

/**
 * @brief Update (compute) the levelset of the grid for a given set of Geometry
 *
 * @param[in] id Unique id of the Geometry
 */
void compute(int idGeo) {
  Geometry * geom = getGeometry(idGeo);
  geom->computeLevelSet(_grid);
}

/**
 * @brief Transport the levelset of the grid for a given set of Geometry
 *
 * @param[in] tag Tag of the Geometrie's set
 * @param[in] u vector
 * @param[in] dt
 */
void transport(std::string tag, const PiercedVector<NPoint> &u, const double dt);

};
}
#endif /* __LEVELSET_HPP__ */
