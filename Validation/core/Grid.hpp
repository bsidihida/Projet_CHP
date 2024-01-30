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

/*
** Neos library
*
* @file Grid.hpp
* @brief This file contains the Grid class
* @author Matias Hastaran
* @version 0.1
* @date 2017-05-14
* @copyright Inria
*/
#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <vector>
#include <string>
#include <cell.hpp>
#include <voloctree.hpp>
#include <adaption.hpp>
#include <common.hpp>
#include <BoundaryConditions.hpp>

namespace neos {

/**
 * @class Grid
 * @brief Abstraction class to use bitpit VolOctree class
 * @details Grid is the main Neos Object, as entrypoint you can set:
 * - Grid size and dimension (2 or 3D)
 * - MPI Communicator
 * - Per face Boundaries Conditions (Neumann, Dirichlet)
 * - And you can choose which cell to refine or to coarse from a given std::vector
 */
class Grid : public bitpit::VolOctree
{
private:
int _dim;                         /**< Grid dimension (2D or 3D) */
BoundaryConditions *_BC;
int _nFaceVertex;
std::vector<std::vector<int> > _faceVertex;
bool haveBC = false;
void setupGrid();
std::vector<bitpit::adaption::Info> _mapper;
/*
Because VTKStreamer use std_unique_ptr, VTK.addData isn't persistent
when using python bindings. So we to cache scalar and array values
before write
*/
std::vector<std::vector<double>> _pybind_vector_cache;
std::vector<std::vector<std::array<double, 3>>> _pybind_array_cache;

public:
/**
 * @brief Grid class constructor, inherit from bitpit::VolOctree constructor
 *
 * @param[in] X X origin
 * @param[in] Y Y origin
 * @param[in] Z Z origin
 * @param[in] L Length of the domain
 * @param[in] dh Maximum allowed cell size of the initial refinement
 * @param[in] dim is the dimension of the grid
 * @param[in] haloSize is the size, expressed in number of layers, of the ghost cells halo
 * @param[in] communicator	is the communicator to be used for exchanging data among the processes. If a null comunicator is provided, a serial patch will be created
 */

  Grid(double X, double Y, double Z, double L, double dh, uint8_t dim = GRID_3D,
     uint8_t haloSize = 1, MPI_Comm communicator = MPI_COMM_WORLD);

/**
 * @brief Grid class constructor from bitpit::PabloUniform
 * @param[in] treePointer is the tree that will be used
 * @param[in] dim is the dimension of the grid
 **/

Grid(std::unique_ptr<bitpit::PabloUniform>& treePointer, uint8_t dim= GRID_3D) :
  bitpit::VolOctree(std::move(treePointer), &treePointer)
{
  _dim = dim;
  haveBC = false;

  setupGrid();
}


/**
 * @brief Overloaded constructor with boundary conditions
 *
 * @param X
 * @param Y
 * @param Z
 * @param L
 * @param dh
 * @param dim
 * @param BC
 *
 * @return
 */
Grid(
  double X, double Y, double Z,
  double L, double dh,
  BoundaryConditions *BC,
  uint8_t dim = GRID_3D,
  uint8_t haloSize = 1,
  MPI_Comm communicator = MPI_COMM_WORLD
  );


/**
  * @brief Grid class destructor
  */
~Grid()
{
}


/**
 * @brief Get the dimension of the Grid
 */
int getDim() {
  return _dim;
}

/**
 * @brief Get the number of Vertices per Face
 */
int getNFaceVertex() {
  return _nFaceVertex;
}

/**
 * @brief gets the number of cells
 * @return Cells numbers (only internal Cells, no ghosts)
 */
long nbCells() {
  return getInternalCellCount();
}


/**
 * @brief Get the smallest edge size
 *
 *
 * @return Smallest edge size
 */
double getMinSize();

/**
 * @brief Refine a Grid from a Vector Range Values
 * @details For a given vector of a Levelset (values in center of each cell),
 * we refine if the value (distance) is between min and max value
 * @param[in] v Vector of Levelset values for each cell
 * @param[in] min Distance value above we refine
 * @param[in] max Distance value below we refine
 * @return The number of refined cells
 */
long refineFromRangeVector(std::vector<double> v, double min, double max);

/**
 * @brief Coarse the Grid from a Vector Range Values
 * @details For a given vector of a Levelset (values in center of each cell),
 * we coarse if the value (distance) is greater than a threshold
 * @param[in] v Vector of Levelset values for each cell
 * @param[in] threshold Distance value  above we coarse the Cell
 * @return The number of coarsed cells
 */
long coarseFromRangeVector(std::vector<double> v, double threshold);

std::vector<bitpit::adaption::Info> getMapper() {
  return _mapper;
}
/**
 * @brief Update Mesh for Refinement and store Mapper for Levelset Updates
 */
 void update(bool trackAdaption = true, bool squeezeStorage = false);

/**
 * @brief Return a vector of the neighbours id
 *
 * @param[i] x0 Coordinate of the point
 * @param[i] id Id of the cell
 */
std::vector<long> findQuadrantNeighs(const NPoint &x0, const long id);

/**
 * @brief Return the id of the cell
 *
 * @param[i] point Coordinate of the point
 */
long locatePoint(const NPoint &point);


/**
 * @brief set name for export files
 *
 * @param[i] name file name
 */
void setExportName(std::string name);

/**
 * @brief start counter for animated VTK file
 *
 * @param[in] c first frame for animation
 */
void setCounter (int c=0) {
  getVTK().setCounter(c);
}

/**
 * @brief reset counter to stop animated VTK file. return last count
 *
 */
int unsetCounter () {
  return getVTK().unsetCounter();
}

/**
 * @brief Add a set of data to export
 *
 * @param[in] tag Tag name of the dataset
 * @param[in] data Vector of data
 */
void addData(std::string tag, std::vector<double> &data)
{
  #ifdef BINDER_METHODS
    _pybind_vector_cache.push_back(data);
    getVTK().addData(tag, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, _pybind_vector_cache.back());
  #else
  // You have to write it (Grid::write)
  getVTK().addData(tag, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL,data);
  #endif
}

/**
 * @brief Add a set of data to export
 *
 * @param[in] tag Tag name of the dataset
 * @param[in] data Vector of array of data
 */
void addData(std::string tag, std::vector<std::array<double, 3> > &data)
{
  #ifdef BINDER_METHODS
    _pybind_array_cache.push_back(data);
    getVTK().addData(tag, bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::CELL, _pybind_array_cache.back());
  #else
  // You have to write it (Grid::write)
  getVTK().addData(tag, bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::CELL, data);
  #endif
}

/**
 * @brief Write VTK file
 *
 */
void write() {
  bitpit::VolOctree::write();
  #ifdef BINDER_METHODS
  _pybind_vector_cache.clear();
  _pybind_array_cache.clear();
  #endif
}

/**
  * @brief Get vertices of a given Face
  * @details For a given face, return vertices connected to this face
  * param[in] Face id
  * 2D face Id numbering: 0=left, 1=right, 2=bottom, 3=up
  * 2D vertex Id numbering: 0=bottom left 1=bottom right 2:upper left 3=upper right
  * 3D face Id numbering: 0=left, 1=right, 2=bottom, 3=up, 4=front, 5=rear
  * 3D vertex Id numbering: 0=front bottom left 1=front bottom right 2=front upper left 3=front upper right
  * 4=rear bottom left 5=rear bottom right 6=rear upper left 7=rear upper right
  * @return List of vertices Ids
  */
std::vector<int> getFaceVertexLocalIds(int face);

bool hasVertexOnBorder(const int& interId,
                       int& interOnBorder);

/**
 * @brief Is the interface "interId" is a domain boundary cell ?
 *
 * @param[in] interId Global id of the interface
 *
 * @return True if the interface is on the domain boundary
 */
bool isInterfaceDomainBorder(const long& interId );

/**
 * @brief Is the interface "inter" a domain boundary cell ?
 *
 * @param[in] inter : the interface
 *
 * @return True if the interface is on the domain boundary
 */
bool isInterfaceDomainBorder(const bitpit::Interface& inter);

/**
 * @brief Is the cell "cellId" is a boundary cell ?
 *
 * @param[in] cellId Global id of the cell
 *
 * @return True if the cell is on boundary
 */
bool isBorder(const long& cellId );

/**
 * @brief Return number of faces which are on boundary
 *
 * @param cellId Global id of the cell
 *
 * @return Number of boundary faces
 */
size_t nbBorder(const long& cellId );

/**
 * @brief Return the local id of the boundary face.
 *
 * @param cellId Global id of the cell
 *
 * @return local id of the boundary face
 */
int getBoundFaceId(const long& cellId);

/**
 * @brief Return the local ids of the boundary faces.
 *
 * @param cellId Global id of the cell
 *
 * @return vector of boundary faces local id
 */
std::vector<int> getBoundFacesId(const long& cellId);

/**
 * @brief compute boundary value depending on boundary conditions (Neumann
 *     or Dirichlet).
 *
 * @param[in] coord: Coordinate where we want to compute boundary value
 * @param[in] val value at cell center
 * @param[in] boundaryFace local id of the boundary face
 * @param[in] type Type of function Var::P, Var::Ux, Var::Uy or Var::LS
 *
 * @return Value used for boundary ghost cell.
 */
double computeBoundaryValue(const NPoint& coord,
                            const double& val,
                            const int& boundaryFace,
                            Var type = Var::LS,
                            double t=0.) const;

/**
 * @brief Return boundary conditions
 *
 *
 * @return BoundaryConditions for the grid
 */
BoundaryConditions* getBoundaryConditions();

/**
 * @brief Does the user provides boundary conditions in constructor ?
 *
 *
 * @return True if user provides boundary conditions. False if not
 */
bool haveBoundaryConditions();

std::unordered_map<long, long> getCellGlobalMap();

#ifdef BINDER_METHODS
/**
  * Those methods a explicitly added for python bindings generation
  *
  */

// PiercedVector<bitpit::Cell>& getCells() {
//   bitpit::PiercedVector<bitpit::Cell> temp = bitpit::VolOctree::getCells();
//   PiercedVector<bitpit::Cell> temp2=temp;
//   return temp2;
// }

const bitpit::PiercedVector<bitpit::Cell>& getCells() const {
  return bitpit::PatchKernel::getCells();
}

const bitpit::Cell& getCell(const long int &id) const {
  return bitpit::VolOctree::getCell(id);
}

int getCellLevel(const long int &id) const {
  return bitpit::VolOctree::getCellLevel(id);
}

// const MPI_Comm & getCommunicator() const
//  {
//      return bitpit::VolOctree::getCommunicator();
//  }

#endif // BINDER_METHODS

};
}


#endif /* _GRID_HPP_ */
