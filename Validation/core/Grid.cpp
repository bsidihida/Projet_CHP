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
 * Neos library
 *
 * @file Grid.cpp
 *
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-11-23
 * @copyright Inria
 */
#if ENABLE_MPI==1
#include <mpi.h>
#endif
#include "Grid.hpp"
#include "Utils.hpp"
#include "NeosPiercedVector.hpp"
#include "NeosAssert.hpp"
namespace neos {

#if ENABLE_MPI==1
  Grid::Grid(double X, double Y, double Z, double L, double dh, uint8_t dim,
    uint8_t haloSize, MPI_Comm communicator)
    : bitpit::VolOctree(-1, dim, { { X, Y, Z} }, L, dh, communicator, haloSize)
{
  _dim = dim;
  haveBC = false;

  setupGrid();
}
#else
Grid::Grid(double X, double Y, double Z, double L, double dh, uint8_t dim, uint8_t haloSize)
  : bitpit::VolOctree(-1, dim, { { X, Y, Z} }, L, dh, haloSize)
{
  _dim = dim;
  haveBC = false;

  setupGrid();
}
#endif

Grid::Grid(
  double X, double Y, double Z,
  double L, double dh,
  BoundaryConditions *BC,
  uint8_t dim,
  uint8_t haloSize,
  MPI_Comm communicator
) : Grid(X,Y,Z,L,dh,dim,haloSize,communicator)
{
  _BC   = BC;
  haveBC = true;
}

// Set a file name and set the counter to 0 (for multiple write)
void Grid::setExportName(std::string name)
{
  getVTK().setName(name);
  setVTKWriteTarget(bitpit::PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
}

long Grid::refineFromRangeVector(std::vector<double> v, double min, double max){
  long counter = 0;

  PiercedVector<double> pv = VtoPV(v, this);

  for (auto &cell : getCells()) {
    const long &id = cell.getId();
    if (getCell(id).isInterior()
        &&
        (pv[id] < max)
        &&
        (pv[id] > min)) {
      markCellForRefinement(id);
      counter++;
      //std::cerr << " raf " << maxr << std::endl;
    }
  }
  update(false,true);
  return counter;
}

long Grid::coarseFromRangeVector(std::vector<double> v, double threshold){
  long counter = 0;

  PiercedVector<double> pv = VtoPV(v, this);

  for (auto &cell : getCells()) {
    const long &id = cell.getId();
    if (getCell(id).isInterior() && (pv[id] > threshold)) {
      markCellForCoarsening(id);
      counter++;
    }
  }
  update(false,true);
  return counter;
}

// Return the neighbours of one point
std::vector<long> Grid::findQuadrantNeighs(const NPoint &x0, const long id)
{
  int vertex = 0;
  NPoint center = evalCellCentroid(id);

  if (x0[NPX] >= center[NPX]) {
    ++vertex;
  }
  if (x0[NPY] >= center[NPY]) {
    vertex += 2;
  }
  if (_dim == GRID_3D && x0[NPZ] >= center[NPZ]) {
    vertex += 4;
  }
  return findCellVertexOneRing(id, vertex);
}

void Grid::update(bool trackAdaption, bool squeezeStorage)
{
  _mapper = bitpit::VolOctree::update(trackAdaption, squeezeStorage);
}

// Abstraction to bitpit, return the cell ID
long Grid::locatePoint(const NPoint &point)
{
  return bitpit::VolOctree::locatePoint((const std::array<double, 3>)point);
}

// // Add new data (in VTK) with a new tag
// void Grid::addData(std::string tag, std::vector<double> &data)
// {
//   #ifdef BINDER_METHODS
//     _vector_cache.push_back(data);
//   #endif
//   // You have to write it (Grid::write)
//   getVTK().addData(tag, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, data);
// }
//
// // Add new data (in VTK) with a new tag
// void Grid::addData(std::string tag, std::vector<std::array<double, 3> > &data)
// {
//   #ifdef BINDER_METHODS
//     _array_cache.push_back(data);
//   #endif
//   // You have to write it (Grid::write)
//   getVTK().addData(tag, bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::CELL, data);
// }
//
// void Grid::write() {
//   bitpit::VolOctree::write();
//   #ifdef BINDER_METHODS
//   _vector_cache.clear();
//   _array_cache.clear();
//   #endif
// }

double Grid::getMinSize()
{
  // Get minimal size on each process
  double dxMin(this->getTree().getLocalMinSize());
  double dxMinLoc(dxMin);

  // If more than on process, get the global minimum.
  if (getProcessorCount()>1)
  {
    MPI_Barrier(this->getCommunicator());
    MPI_Allreduce(&dxMinLoc,
                  &dxMin,
                  1,
                  MPI_DOUBLE,
                  MPI_MIN,
                  this->getCommunicator());
  }
  return dxMin;
}

std::vector<int> Grid::getFaceVertexLocalIds(int face)
{
  return _faceVertex[face];
}

bool Grid::hasVertexOnBorder(const int& interId,
                             int& interOnBorder)
{
  bitpit::Interface& inter = this->getInterface(interId);
  bitpit::Interface crossInter;
  bitpit::Cell& cell = this->getCell(inter.getOwner());

  const int& cellId = cell.getId();
  int face = inter.getOwnerFace();
  bool hasVertexOnBord(false);

  if (inter.isBorder())
  {
    interOnBorder = face;
    return true;
  }
  else if (this->isBorder(cellId))
  {
    std::vector<int> boundFaces = this->getBoundFacesId(cellId);

    switch(face/2)
    {
    case 0:
      if (std::find(boundFaces.begin(), boundFaces.end(), 2) !=
          boundFaces.end())
      {
        hasVertexOnBord = true;
        interOnBorder   = 2;
      }

      if (std::find(boundFaces.begin(), boundFaces.end(), 3) !=
          boundFaces.end())
      {
        hasVertexOnBord = true;
        interOnBorder   = 3;
      }
      break;
    case 1:
      if (std::find(boundFaces.begin(), boundFaces.end(), 0)
          != boundFaces.end())
      {
        hasVertexOnBord = true;
        interOnBorder   = 0;
      }

      if (std::find(boundFaces.begin(), boundFaces.end(), 1) !=
          boundFaces.end())
      {
        hasVertexOnBord = true;
        interOnBorder   = 1;
      }
      break;
    default:
      break;
    }
  }

  return hasVertexOnBord;
}

bool Grid::isInterfaceDomainBorder(const long& id)
{
  return isInterfaceDomainBorder(this->getInterface(id));
}

bool Grid::isInterfaceDomainBorder(const bitpit::Interface& inter)
{
  bitpit::Cell& cell = this->getCell(inter.getOwner());

  if (inter.isBorder() && cell.isInterior())
  {
    return true;
  }
  return false;
}

bool Grid::isBorder(const long& id)
{

  for (int i=0; i<4; i++)
  {
    if (this->findCellFaceNeighs(id, i).size() == 0)
      return true;
  }
  return false;
}

size_t Grid::nbBorder(const long& id)
{
  const bitpit::Cell cell = this->getCell(id);
  size_t nbBound(0);
  if (cell.isInterior())
  {
    for (size_t i = 0; i < 4; i++)
    {
      if (this->findCellFaceNeighs(id, i).size() == 0)
        nbBound++;
    }
  }
  return nbBound;
}

int Grid::getBoundFaceId(const long& id)
{
  const bitpit::Cell cell = this->getCell(id);
  if (!this->isBorder(id))
  {
    std::cerr<<"This is not a border cell"<<std::endl;
    return -1.;
  }
  else
  {
    for (size_t i = 0; i < 4; i++)
    {
      if (cell.isFaceBorder(i))
        return i;
    }
  }
  ASSERT(false, "Unreachable statement");
  return -1;
}

std::vector<int> Grid::getBoundFacesId(const long& id)
{
  const bitpit::Cell cell = this->getCell(id);
  std::vector<int> boundFaces;
  if (!this->isBorder(id))
  {
    std::cerr<<"This is not a border cell"<<std::endl;
  }
  else
  {
    for (size_t i = 0; i < 4; i++)
    {
      if (cell.isFaceBorder(i))
        boundFaces.push_back(i);
    }
  }
  return boundFaces;
}

//double Grid::computeBoundaryValue(const int &cellId, const double& val, const int& boundaryFace)
//{
//    std::unordered_map<int, std::pair<std::string, std::array<double,3>>>
//	cond = _BC->getBCConditions();
//
//    if ( cond[boundaryFace].first == "Dirichlet" )
//    {
//
//	// Extrapol value on ghost boundary cell
//	// FIXME: This condition works only on rectangular domains.
//	return 2*cond[boundaryFace].second[boundaryFace/2] - val;
//    }
//    else if (cond[boundaryFace].first == "Neumann")
//    {
//	return val;
//    }
//}

double Grid::computeBoundaryValue(const NPoint& coord,
                                  const double& val,
                                  const int& boundaryFace,
                                  Var type,
                                  double t) const
{

  std::pair<Var, int> varFace = std::make_pair(type, boundaryFace);
  auto& cond = _BC->getBCConditions();

  if ( cond[varFace].first == "Dirichlet" )
  {
    return cond[varFace].second(coord, t);
  }
  else if (cond[varFace].first == "Neumann")
  {
    return val;
  }
  ASSERT(false, "Unreachable statement");
  return -1.0;
}

BoundaryConditions* Grid::getBoundaryConditions()
{
  return _BC;
}

bool Grid::haveBoundaryConditions()
{
  return haveBC;
}

std::unordered_map<long, long> Grid::getCellGlobalMap()
{
  return bitpit::PatchNumberingInfo(this).getCellGlobalMap();
}

void Grid::setupGrid()
{
  initializeAdjacencies();
  initializeInterfaces();
  update();
#if defined (ENABLE_MPI)
  if (getProcessorCount()>1 && isCommunicatorSet()) {
    partition(true);
  }
#endif

  int face = 2 * _dim;
  _nFaceVertex = face - 2;
  _faceVertex.resize(face, std::vector<int>(_nFaceVertex, 0));

  /* 2D:
   * segment: 0=gauche, 1=droite, 2=bas, 3=haut
   * points:  0=bas_gauche 1=bas_droite 2:haut_gauche 3=haut_droite
   * facevertex[ numero_face ] = liste des sommets de ce segment
   *
   * 3D : idem avec des faces et les sommets des faces
   */

  if (_dim == 2)
  {
    _faceVertex[0] = {{ 0, 2 }};
    _faceVertex[1] = {{ 1, 3 }};
    _faceVertex[2] = {{ 0, 1 }};
    _faceVertex[3] = {{ 2, 3 }};
  }
  else
  {
    _faceVertex[0] = {{ 0, 2, 4, 6 }};
    _faceVertex[1] = {{ 1, 3, 5, 7 }};
    _faceVertex[2] = {{ 0, 1, 4, 5 }};
    _faceVertex[3] = {{ 2, 3, 6, 7 }};
    _faceVertex[4] = {{ 0, 1, 2, 3 }};
    _faceVertex[5] = {{ 4, 5, 6, 7 }};
  }
}
}
