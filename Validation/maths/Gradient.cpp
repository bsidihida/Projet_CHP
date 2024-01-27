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
 * @file Gradient.cpp
 *
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 * @copyright Inria
 */

#include <cell.hpp>
#include <interface.hpp>
#include "Gradient.hpp"
#include "InterpolatorFactory.hpp"
#include "Utils.hpp"
#include "Function.hpp"
#include "MLS2D.hpp"
#include "NeosAssert.hpp"
#include <numeric>

namespace neos {

UserDataComm<double>* Gradient::_userComm;
int Gradient::_userComm_only_one_init_dirty=0;

Gradient::Gradient(int dim, Grid *grid, StencilBuilder* stencil)
  : _grid(grid),
  _stencil(stencil),
  _hasStencil(stencil != nullptr)
{
  _dim = dim;
  _mapGlobal = _grid->getCellGlobalMap();
}

/* param:
 * OwnerId: OwnerId d'une cellule,
 * localVertex: l'indice d'un sommet
 * */

std::vector<std::pair<long,double> > Gradient::buildVertexStencil(const long ownerId, const int localVertex) {

  NPoint x0     = _grid->getVertexCoords(_grid->getCell(ownerId).getVertexId(localVertex));                // vertex coordinate
  std::vector<long>   neighbours = _grid->findCellVertexOneRing(ownerId, localVertex);   // cell list around this vertex
  std::vector<std::array<double,3> > xref;


  for (std::size_t iNeigh = 0; iNeigh<neighbours.size(); ++iNeigh)   // for each neibhour cell id
    xref.push_back(_grid->evalCellCentroid(neighbours[iNeigh]));     // keep the centroid

  IInterpolator *interp = InterpolatorFactory::get(_iType);
  std::vector<double> weights = interp->computeInterpolation(x0, xref);   // get weights by interpolation

  std::vector<std::pair<long,double> > stencil;
  for (std::size_t i=0; i<neighbours.size(); ++i) { // for each neibhour cell id
    std::pair<long,double> temp(neighbours[i],weights[i]);     // associate <cell id, weights>
    stencil.push_back(temp);
  }
  return stencil;   // return a map <cell id, weight>
}



std::vector<std::pair<long,double> > Gradient::computeNormalGradient(const long interfaceId, std::vector<std::vector<std::pair<long,double> > > &grad){

  std::array<long,2> owners = _grid->getInterface(interfaceId).getOwnerNeigh();
  std::array<double,3> normal = _grid->evalInterfaceNormal(interfaceId);
  int dir = fabs(normal[1])+2*fabs(normal[2]);

  std::array<double,3> centerVec, c1 = _grid->evalCellCentroid(owners[1]), c0 = _grid->evalCellCentroid(owners[0]);
  double d=0;
  for (std::size_t i=0; i<3; ++i) {
    centerVec[i] = c1[i]-c0[i];
    d += (c1[i]-c0[i])*(c1[i]-c0[i]);
  }
  double normInv = 1/sqrt(d);

  for (std::size_t i = 0; i<3; ++i)
    centerVec[i] *= normInv;

  d = 1.0 / DotProduct(centerVec, normal);

  std::vector<std::pair<long,double> > normalGradient;
  normalGradient.push_back(std::make_pair(owners[0],-d*normInv));
  normalGradient.push_back(std::make_pair(owners[1],d*normInv));


  if (_grid->evalCellVolume(owners[0]) != _grid->evalCellVolume(owners[1])) {
    int component = 0;
    if (grad.size() !=0 ) {
      for (int i=0; i<_dim; ++i) {
        if (i!=dir) {

          for ( std::size_t j = 0; j < grad[component].size(); ++j) {

            grad[component][j].second = -grad[component][j].second*centerVec[i] * d;

          }

          ++component;

        }
      }
    }

    for (std::size_t i = 0; i < grad.size(); ++i) {

      for (std::size_t j = 0; j < grad[i].size(); ++j) {

        bool exist = false;
        for (std::size_t k = 0; k < normalGradient.size(); ++k) {

          if (normalGradient[k].first == grad[i][j].first) {
            normalGradient[k].second += grad[i][j].second;
            exist = true;
          }

        }
        if (!exist) normalGradient.push_back(grad[i][j]);

      }

    }
  }

  //for (std::size_t i = 0; i < normalGradient.size();++i)
  //	normalGradient[i].second *= area;

  return normalGradient;

}

// for a given interface, return a map <vertex id, associated "size">
std::map<std::pair<int,int>,double> Gradient::buildSizeMap(const long interfaceId){

  long ownerFace = _grid->getInterface(interfaceId).getOwnerFace();   // face index (3D: from 0 to 6)
  long ownerId = _grid->getInterface(interfaceId).getOwner();   // owner cell Id

  double interfaceSizeInv = 1;
  if (_dim == 2)
    interfaceSizeInv = 1/_grid->evalInterfaceArea(interfaceId);     // inverse of length in 2D
  else if (_dim == 3)
    interfaceSizeInv = 1/sqrt(_grid->evalInterfaceArea(interfaceId));     // inverse of surface in 3D

  std::map<std::pair<int,int>,double> sizeMap;

  bool isbound = false;

  for (int vertex = 0; vertex < _grid->getNFaceVertex(); ++vertex) {

    std::vector<std::pair<long,double> > vertexStencil = buildVertexStencil(ownerId,_grid->getFaceVertexLocalIds(ownerFace)[vertex]);    // weights unused ?? TODO

    if (vertexStencil.size() == std::pow(2, _dim-1)) {

      isbound = true;
    }

    sizeMap[{vertex,0}] = interfaceSizeInv;

  }

  if (isbound) {
    for (int vertex = 0; vertex < _grid->getNFaceVertex(); ++vertex) {
      sizeMap[{vertex,0}] = interfaceSizeInv * 2;
    }
  }

  return sizeMap;

}



std::vector<std::vector<std::pair<long,double> > > Gradient::buildInterfaceGradientStencil(const long interfaceId){

  long ownerFace = _grid->getInterface(interfaceId).getOwnerFace();
  long ownerId = _grid->getInterface(interfaceId).getOwner();
  bitpit::Cell ownerCell = _grid->getCell(ownerId);

  double interfaceSizeInv = 1;
  if (_dim == 2)
    interfaceSizeInv = 1/_grid->evalInterfaceArea(interfaceId);
  else if (_dim == 3)
    interfaceSizeInv = 1/sqrt(_grid->evalInterfaceArea(interfaceId));

  std::map<std::pair<int,int>,double> invSize = buildSizeMap(interfaceId);

  std::vector<std::vector<std::pair<long,double> > > grad;

  for (int vertex = 0; vertex < _grid->getNFaceVertex(); ++vertex) {  // for each interface vertex 2D: 0-1  3D: 0-3

    std::vector<std::pair<long,double> > vertexStencil = buildVertexStencil(ownerId, _grid->getFaceVertexLocalIds(ownerFace)[vertex]);


    if (_dim==2) {
      if ((vertex+1)%2==0) {      // if right side vertex (odd index)
        int dir = vertex/2;         // dir=0:bottom vertex    dir=1:top vertex
        for (std::size_t ivS = 0; ivS < vertexStencil.size(); ++ivS)
        {
          bool exist = false;
          for (std::size_t igrad=0; igrad < grad[dir].size(); ++igrad)
          {

            if (vertexStencil[ivS].first == grad[dir][igrad].first)
            {
              grad[dir][igrad].second += vertexStencil[ivS].second * invSize[{vertex,0}];
              exist = true;
            }

          }
          if (!exist)
          {
            grad[dir].push_back(vertexStencil[ivS]);
            grad[dir][grad[dir].size()-1].second *= invSize[{vertex,0}];
          }
        }
      }
      else{       // if left side vertex (even index)
        for (std::size_t i=0; i < vertexStencil.size(); ++i)
          vertexStencil[i].second *= -invSize[{vertex,0}];
        grad.push_back(vertexStencil);
      }
    }

    else if (_dim==3) {

      for (std::size_t ivS = 0; ivS < vertexStencil.size(); ++ivS) {

        for (unsigned int dir = 0; dir < 2; ++dir) {
          bool exist = false;
          if (grad.size() > dir ) {
            for (std::size_t igrad=0; igrad < grad[dir].size(); ++igrad) {

              if (vertexStencil[ivS].first == grad[dir][igrad].first) {

                if ( ( dir == 0 && (vertex == 0 || vertex == 2)) ||
                     ( dir == 1 && (vertex == 0 || vertex == 1)) ) grad[dir][igrad].second -= vertexStencil[ivS].second*interfaceSizeInv;

                if ( ( dir == 0 && (vertex == 1 || vertex == 3)) ||
                     ( dir == 1 && (vertex == 2 || vertex == 3)) ) grad[dir][igrad].second += vertexStencil[ivS].second*interfaceSizeInv;
                exist = true;
              }

            }
          }

          if (!exist) {
            if (grad.size() > dir) {
              grad[dir].push_back(vertexStencil[ivS]);
              if ( ( dir == 0 && (vertex == 0 || vertex == 2)) ||
                   ( dir == 1 && (vertex == 0 || vertex == 1)) ) grad[dir][grad[dir].size()-1].second *= -interfaceSizeInv;

              if ( ( dir == 0 && (vertex == 1 || vertex == 3)) ||
                   ( dir == 1 && (vertex == 2 || vertex == 3)) ) grad[dir][grad[dir].size()-1].second *= interfaceSizeInv;
            }
            else{
              grad.push_back(std::vector<std::pair<long,double> >(1,vertexStencil[ivS]));
              if ( ( dir == 0 && (vertex == 0 || vertex == 2)) ||
                   ( dir == 1 && (vertex == 0 || vertex == 1)) ) grad[dir][0].second *= -interfaceSizeInv;

              if ( ( dir == 0 && (vertex == 1 || vertex == 3)) ||
                   ( dir == 1 && (vertex == 2 || vertex == 3)) ) grad[dir][0].second *= interfaceSizeInv;
            }

          }
        }
      }

    }

  }

  //in 3D the gradient in the interface is computed on each edge. Take the mean.
  if (_dim == 3) {
    for (unsigned int i = 0; i < grad.size(); ++i) {
      for (unsigned int j = 0; j < grad[i].size(); ++j) {
        grad[i][j].second *= 0.5;
      }
    }
  }

  return grad;
}


std::vector<std::pair<long,double> > Gradient::buildInterfaceNormalGradientStencil(const long interfaceId){

  std::vector<std::pair<long,double> > m_normalGradientWeights;
  std::vector<std::pair<long,double> > m_normalGradientInterfaceId;

  // if (m_normalGradientWeights.size() == 0 || m_normalGradientInterfaceId != interfaceId ){
  //  m_normalGradientInterfaceId = interfaceId;
  // m_normalGradientWeights.clear();

  //long ownerFace = grid->getInterface(interfaceId).getOwnerFace();
  long ownerId = _grid->getInterface(interfaceId).getOwner();
  long neighId = _grid->getInterface(interfaceId).getNeigh();
  bitpit::Cell ownerCell = _grid->getCell(ownerId);

  double interfaceSizeInv = 1;
  if (_dim == 2)
    interfaceSizeInv = 1/_grid->evalInterfaceArea(interfaceId);
  else if (_dim == 3)
    interfaceSizeInv = 1/sqrt(_grid->evalInterfaceArea(interfaceId));
  std::vector<std::vector<std::pair<long,double> > > grad;
  if (_grid->evalCellVolume(ownerId) == _grid->evalCellVolume(neighId)) {
    //std::vector<std::pair<long,double>> ngrad;
    //ngrad.push_back(std::make_pair(ownerId, -interfaceSizeInv));
    //ngrad.push_back(std::make_pair(neighId, interfaceSizeInv));
    //return ngrad;
    m_normalGradientWeights.push_back(std::make_pair(ownerId, -interfaceSizeInv));
    m_normalGradientWeights.push_back(std::make_pair(neighId,  interfaceSizeInv));
  }
  else{
    grad = buildInterfaceGradientStencil(interfaceId);
    //std::vector<std::pair<long,double>> ngrad = computeNormalGradient(grid, interfaceId, grad);
    //return ngrad;
    m_normalGradientWeights = computeNormalGradient(interfaceId, grad);
  }
  //}
  //else if (m_interfaceId != interfaceId){
  //	throw CRAFPACK::common::GeneralException("Try to use existing normal gradient stencil on wrong interface");
  //}

  return m_normalGradientWeights;
}

std::vector<std::array<double,3> >  Gradient::computeFVGradient(PiercedVector<double> &val)
{
  PiercedVector<std::array<double,3> > gradient;

  for (auto &cell : _grid->getCells()) {
    //if (cell.isInterior()) {
    const long &i = cell.getId();
    gradient.emplace(i);
    gradient[i] = {{0, 0, 0}};
    //}
  }

  for (auto &inter : _grid->getInterfaces())
  {
    const long &interfaceId = inter.getId();
    if (!inter.isBorder())
    {
      std::array<long, 2> owners = inter.getOwnerNeigh();
      std::vector<std::pair<long,double> > weights;
      double area = _grid->evalInterfaceArea(interfaceId);

      weights = buildInterfaceNormalGradientStencil(interfaceId);

      for (std::size_t i = 0; i<weights.size(); ++i)
      {
        const long id = weights[i].first;

        for (int j = 0; j<_grid->getDimension(); ++j)
        {
          gradient[owners[0]][j] += val[id]*area*weights[i].second*_grid->evalInterfaceNormal(interfaceId)[j]/_grid->evalCellVolume(owners[0]);


          gradient[owners[1]][j]-= val[id]*area*weights[i].second*_grid->evalInterfaceNormal(interfaceId)[j]/_grid->evalCellVolume(owners[1]);
        }
      }
    }
  }

  return PVtoV(gradient,_grid);
}

std::vector<std::array<double,3> > Gradient::computeLSGradient(PiercedVector<double> &val){
  PiercedVector<std::array<double,3> > gradient;

  if (Gradient::_userComm_only_one_init_dirty!=10101)
  {
    Gradient::_userComm = new UserDataComm<double>(*_grid, 1);
    Gradient::_userComm_only_one_init_dirty=10101;
  }

  Gradient::_userComm->update(val.size());
  Gradient::_userComm->communicate(val);

  for (auto &cell : _grid->getCells()) {

    const long &id=cell.getId();
    gradient.emplace(id);

    if (_grid->getCell(id).isInterior()) {
      gradient[id]=computeLSCellGradient(id,val);
    }
  }
  return PVtoV(gradient,_grid);
}

std::vector<std::array<double,3> > Gradient::computeLSGradient(const std::vector<double> &val){

  PiercedVector<double> a = VtoPV(val,_grid);
  return computeLSGradient(a);
}

std::array<double,3> Gradient::computeLSCellGradient(long id, PiercedVector<double> &data){

  int spacedim=_grid->getDimension();
  std::vector<long> neigh=_grid->findCellNeighs(id);

  double A[spacedim*spacedim];
  double F[spacedim];

  std::fill_n(A,spacedim*spacedim,0.0);
  std::fill_n(F,spacedim,0.0);

  buildLSCellGridMatrix(id,neigh,A);
  buildLSCellRHS(id,neigh,data,F);

  inverse(spacedim,A);

  std::array<double,3> grad={{0.0,0.0,0.0}};

  for (int i=0; i<spacedim; ++i) {
    for (int j=0; j<spacedim; ++j) {
      grad[i]+=A[i*spacedim+j]*F[j];
    }
  }
  return grad;
}

void Gradient::buildLSCellGridMatrix(long id, std::vector<long> &neigh,double* A){

  int spacedim=_grid->getDimension();
  std::array<double,3> c_id=_grid->evalCellCentroid(id), c_n;

  for (uint32_t i=0; i<neigh.size(); ++i) {
    c_n=_grid->evalCellCentroid(neigh[i]);
    std::array<double,3> dc=c_n-c_id;
    double nd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];
    nd=1.0/nd;

    for (int j=0; j<spacedim; ++j) {
      for (int k=0; k<spacedim; ++k) {
        A[j*spacedim+k]+=nd*dc[j]*dc[k];
      }
    }
  }
}

void Gradient::buildLSCellRHS(long id, std::vector<long> &neigh,PiercedVector<double> &data, double *F){

  int spacedim=_grid->getDimension();
  std::array<double,3> c_id=_grid->evalCellCentroid(id), c_n;

  for (uint32_t i=0; i<neigh.size(); ++i) {
    c_n=_grid->evalCellCentroid(neigh[i]);
    std::array<double,3> dc=c_n-c_id;
    double nd=dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2];
    nd=1.0/nd;

    double df=data[neigh[i]]-data[id];

    for (int j=0; j<spacedim; ++j)
      F[j] += nd*df*dc[j];
  }

}

PiercedVector<std::array<double, 3> > Gradient::computeCCLSGradient(
  PiercedVector<double>& val,
  const double t,
  const Var type)
{
  // FIXME: It may be better to define a template method here.
  auto valVector = std::vector<PiercedVector<double> >{val};
  auto gradients = this->computeCCLSGradient(valVector,
                                             std::vector<Var>{type},
                                             std::vector<double>{t});

  return gradients[0];
}

std::vector<PiercedVector<std::array<double, 3> > > Gradient::computeCCLSGradient(
  std::vector<PiercedVector<double> >& val,
  const std::vector<Var>& types,
  const std::vector<double>& times)
{

  // Define useful variables
  size_t nbVal = val.size();
  std::vector<PiercedVector<std::array<double,3> > > gradients(nbVal);

  // Define communicator if it doesn't exist
  if (Gradient::_userComm_only_one_init_dirty!=10101)
  {
    Gradient::_userComm = new UserDataComm<double>(*_grid, 1);
    Gradient::_userComm_only_one_init_dirty=10101;
  }

  // Communicate ghost values
  for (size_t i=0; i<nbVal; i++)
  {
    Gradient::_userComm->update();
    Gradient::_userComm->communicate(val[i]);
  }

  // Get stencil by reference
  // Careful. If the stencil was not created before, stencils will be a
  // nullptr
  auto& stencils = (_stencil->getCellGradientStencil());



  Stencil sten;

  for (auto &cell: _grid->getCells())
  {
    const long& cellId = cell.getId();

    //If stencil was created, get it. If not build it.
    if (_hasStencil)
      sten = stencils[cellId];
    else
      sten = this->buildCellGradientStencil(cellId);

    for (std::size_t k = 0; k < nbVal; k++)
    {
      gradients[k].emplace(cellId);
      gradients[k][cellId] = {0., 0., 0.};
    }

    // Construct gradient for interior cells.
    // May be the gradient could be construct also for ghost cell, and we
    // could avoid communication at the end ?
    if (_grid->getCell(cellId).isInterior())
    {
      std::size_t nNeighs = sten.neighs.size();
      for (std::size_t i = 0; i < nNeighs; i++)
      {
        const long id = sten.neighs[i];
        std::vector<double> value(nbVal);
        if(sten.isOutsideCell[i] < 0)         // Ghost boundary cell
        {
          int boundaryFace = std::abs(sten.isOutsideCell[i] + 1);

          for (std::size_t k = 0; k < nbVal; k++)
          {
            double t = times[k];
            if (types[k] != Var::LS)
            {
              NPoint ghostOctCenter = _grid->evalCellCentroid(id);
              ghostOctCenter[boundaryFace/2] += pow(-1,boundaryFace+1)
                                                * _grid->evalCellSize(id);

              value[k] = _grid->computeBoundaryValue(ghostOctCenter,
                                                     val[k][id],
                                                     boundaryFace,
                                                     types[k],
                                                     t);
            }
            else
              value[k] = val[k][id] + _grid->evalCellSize(id);
          }
        }
        else
        {
          for (std::size_t k = 0; k<nbVal; k++)
            value[k] = val[k][id];
        }

        // grad[cellId] = sum_i w_i val[neighId]
        for (std::size_t k = 0; k < nbVal; k++)
        {
          gradients[k][cellId][0] += value[k] * sten.weights.weights_dx[i];
          gradients[k][cellId][1] += value[k] * sten.weights.weights_dy[i];
        }
      }
    }
  }

  for (size_t i=0; i<nbVal; i++)
  {
    Gradient::_userComm->update(3);
    Gradient::_userComm->communicate(gradients[i]);
  }
  return gradients;
}

Stencil Gradient::buildCellGradientStencil(const long& cellId)
{
  bool isUniform(true);
  int cellLevel = _grid->getCellLevel(cellId);
  Stencil stencil;
  NPoint cellCenter = _grid->evalCellCentroid(cellId);
  MLS2D interp;
  size_t nbBound = _grid->nbBorder(cellId);

  // Get all neighbours of the cell
  auto neighs = _grid->findCellNeighs(cellId);
  for (auto& nId: neighs)
  {
    if (cellLevel != _grid->getCellLevel(nId))
    {
      isUniform = false;
      break;
    }
  }

  // Non border cell
  if (!_grid->isBorder(cellId))
  {
    //If uniform, we construct weights and neighbours using centered finite
    // difference.
    if (isUniform)
    {
      double h = _grid->evalCellSize(cellId);

      for (int i=0; i< 2*_dim; i++)
      {
        long neighId = _grid->findCellFaceNeighs(cellId, i)[0];
        stencil.neighs.push_back(neighId);
        stencil.isOutsideCell.push_back(0);

        switch(i/2)
        {
        case 0:
          stencil.weights.weights_dx.push_back(pow(-1,i+1) * 0.5/h);
          stencil.weights.weights_dy.push_back(0);
          break;
        case 1:
          stencil.weights.weights_dx.push_back(0);
          stencil.weights.weights_dy.push_back(pow(-1,i+1) * 0.5/h);
          break;
        default:
          break;
        }
      }
    }
    else
    {
      // If non uniform, use a Taylor expansion leading to Least-Square
      // interpolation .
      neighs.insert(neighs.begin(),cellId);
      std::vector<NPoint> xref;
      size_t i(0);
      for (auto &nId: neighs)
      {
        xref.push_back(_grid->evalCellCentroid(nId));
        i++;
      }

      // Compute weights of interpolation
      std::vector<double> weights =
        interp.computeInterpolation(cellCenter,
                                    xref);
      size_t nNeighs = neighs.size();
      for (size_t i=0; i<nNeighs; i++)
      {
        stencil.neighs.push_back(neighs.at(i));
        stencil.isOutsideCell.push_back(0);

        stencil.weights.weights_dx.push_back(weights.at(i * 6 + 1));
        stencil.weights.weights_dy.push_back(weights.at(i * 6 + 2));
      }
    }
  }
  else if (nbBound != 2)   // Border cell with less than 2 face on border.
  {
    if (isUniform)
    {
      // Uniform: centered finite difference
      double h = _grid->evalCellSize(cellId);
      for (int i=0; i<2 * _dim; i++)
      {
        std::vector<long> neighId = _grid->findCellFaceNeighs(cellId,i);

        if (neighId.size() != 0)
        {
          stencil.neighs.push_back(neighId[0]);
          stencil.isOutsideCell.push_back(0);

          switch(i/2)
          {
          case 0:
            stencil.weights.weights_dx.push_back(pow(-1,i+1)*0.5/h);
            stencil.weights.weights_dy.push_back(0);
            break;
          case 1:
            stencil.weights.weights_dx.push_back(0);
            stencil.weights.weights_dy.push_back(pow(-1, i+1) * 0.5/h);
            break;
          default:
            break;
          }
        }
        else
        {
          stencil.neighs.push_back(cellId);
          stencil.isOutsideCell.push_back(-i-1);

          switch(i/2)
          {
          case 0:
            stencil.weights.weights_dx.push_back(pow(-1,i+1)*0.5/h);
            stencil.weights.weights_dy.push_back(0);
            break;
          case 1:
            stencil.weights.weights_dx.push_back(0);
            stencil.weights.weights_dy.push_back(pow(-1,i+1)*0.5/h);
            break;
          default:
            break;
          }
        }
      }
    }
    else
    {
      // Non uniform: Least-square interpolation
      neighs.insert(neighs.begin(),cellId);
      std::vector<NPoint> xref;
      std::vector<int> isGhostCell;
      for (auto &nId: neighs)
      {
        xref.push_back(_grid->evalCellCentroid(nId));
        isGhostCell.push_back(0);
      }

      // Add now the three outside ghost cells for boundary contribution
      const int borderFace = _grid->getBoundFaceId(cellId);
      long vertexNeigh;
      NPoint xBuff(cellCenter);

      xBuff[borderFace/2] += pow(-1, borderFace + 1)
                             * _grid->evalCellSize(cellId);

      xref.push_back(xBuff);
      neighs.push_back(cellId);
      isGhostCell.push_back(-borderFace - 1);

      for (auto &vertex: _grid->getFaceVertexLocalIds(borderFace))
      {
        vertexNeigh = _grid->findCellVertexNeighs(cellId, vertex)[0];

        xBuff = _grid->evalCellCentroid(vertexNeigh);
        xBuff[borderFace/2] += pow(-1, borderFace + 1)
                               * _grid->evalCellSize(vertexNeigh);


        xref.push_back(xBuff);

        neighs.push_back(vertexNeigh);
        isGhostCell.push_back(-borderFace - 1);
      }

      // Compute weights
      std::vector<double> weights =
        interp.computeInterpolation(cellCenter,
                                    xref);

      size_t nNeighs = neighs.size();
      for (size_t i=0; i<nNeighs; i++)
      {
        stencil.neighs.push_back(neighs.at(i));
        stencil.isOutsideCell.push_back(isGhostCell.at(i));

        stencil.weights.weights_dx.push_back( weights.at(i * 6 + 1));
        stencil.weights.weights_dy.push_back( weights.at(i * 6 + 2));
      }
    }
  }
  else   // Corner cell: Least-square interpolation
  {
    neighs.insert(neighs.begin(), cellId);
    std::vector<NPoint> xref;
    std::vector<int> isGhostCell;
    for (auto &nId: neighs)
    {
      xref.push_back(_grid->evalCellCentroid(nId));
      isGhostCell.push_back(0);
    }

    // Add now the four cells for boundary contribution
    const std::vector<int> borderFaces = _grid->getBoundFacesId(cellId);
    for (auto &borderF: borderFaces)
    {
      std::vector<long> vertexNeigh;

      NPoint xBuff(cellCenter);

      xBuff[borderF/2] += pow(-1, borderF + 1)
                          * _grid->evalCellSize(cellId);

      xref.push_back(xBuff);
      neighs.push_back(cellId);
      isGhostCell.push_back(-borderF - 1);
      for (auto &vertex: _grid->getFaceVertexLocalIds(borderF))
      {
        vertexNeigh = _grid->findCellVertexNeighs(cellId, vertex);
        if (vertexNeigh.size() != 0)
        {
          xBuff = _grid->evalCellCentroid(vertexNeigh[0]);
          xBuff[borderF/2] += pow(-1, borderF + 1)
                              * _grid->evalCellSize(vertexNeigh[0]);

          xref.push_back(xBuff);
          neighs.push_back(vertexNeigh[0]);
          isGhostCell.push_back(-borderF - 1);
        }
      }
    }

    std::vector<double> weights =
      interp.computeInterpolation(cellCenter,
                                  xref);

    size_t nNeighs = neighs.size();

    for (size_t i=0; i<nNeighs; i++)
    {
      stencil.neighs.push_back(neighs.at(i));
      stencil.isOutsideCell.push_back(isGhostCell.at(i));

      stencil.weights.weights_dx.push_back( weights.at(i * 6 + 1));
      stencil.weights.weights_dy.push_back( weights.at(i * 6 + 2));

    }
  }
  return stencil;
}

PiercedVector<double> Gradient::computeFCLSGradient(
  PiercedVector<double>& val,
  const double t,
  const Var type)
{
  // FIXME: It may be better to define a template method here.

  auto valVector = std::vector<PiercedVector<double> >{val};

  auto gradients = this->computeFCLSGradient(valVector,
                                             std::vector<Var>{type},
                                             std::vector<double>{t});

  //PiercedVector<double> gradient;
  //for (auto &id: gradients[0].getIds())
  //{
  //	gradient.emplace(id);
  //	gradient[id] = gradients[0][id];
  //}

  return gradients[0];
}

std::vector<PiercedVector<double> > Gradient::computeFCLSGradient(
  std::vector<PiercedVector<double> >& val,
  const std::vector<Var>& types,
  const std::vector<double>& times)
{
  // FIXME: has to be encapsulate in a debug macro
  ASSERT(types.size() == val.size(),
         "Vector of types for boundary conditions has not the same size as \
    vector of values in Gradient::computeFCLSGradient method. You have to \
    specify a type of boundary conditions for each vector you want to compute \
    gradient.");

  size_t nbVal = val.size();

  std::vector<PiercedVector<double> > gradients(nbVal);
  //fcGradientStencil stencil;
  Stencil stencil;

  if (Gradient::_userComm_only_one_init_dirty!=10101)
  {
    Gradient::_userComm = new UserDataComm<double>(*_grid, 1);
    Gradient::_userComm_only_one_init_dirty=10101;
  }

  for (size_t i=0; i<nbVal; i++)
  {
    Gradient::_userComm->update();
    Gradient::_userComm->communicate(val[i]);
  }

  auto& stencils = (_stencil->getInterfaceGradientStencil());
  Stencil sten;

  for (auto &inter: _grid->getInterfaces())
  {
    const long& interId = inter.getId();

    if (_hasStencil)
      sten = stencils[interId];
    else
      sten = this->buildFaceNormalGradientStencil(interId);

    for(size_t i=0; i<nbVal; i++)
    {
      gradients[i].emplace(interId);
      gradients[i][interId] = 0.;
    }


    std::size_t nNeighs = sten.neighs.size();

    for (std::size_t i = 0; i < nNeighs; i++)
    {
      const long id = sten.neighs[i];
      std::vector<double> value(nbVal);

      if(sten.isOutsideCell[i] < 0)       //Ghost boundary cell
      {
        int boundaryFace = std::abs(sten.isOutsideCell[i] + 1 );
        for (std::size_t k = 0; k < nbVal; k++)
        {
          double t = times[k];
          if (types[k] != Var::LS)
          {

            NPoint ghostOctCenter =  _grid->evalCellCentroid(id);
            ghostOctCenter[boundaryFace/2] += pow(-1,boundaryFace+1)
                                              * _grid->evalCellSize(id);

            value[k] = _grid->computeBoundaryValue(ghostOctCenter,
                                                   val[k][id],
                                                   boundaryFace,
                                                   types[k],
                                                   t);
          }
          else
          {
            value[k] = val[k][id] + _grid->evalCellSize(id);
          }
        }
      }
      else
      {
        for (std::size_t k = 0; k < nbVal; k++)
        {
          value[k] = val[k][id];
        }
      }
      for (std::size_t k = 0; k < nbVal; k++)
      {
        gradients[k][interId] += value[k] *
                                 sten.weights.weights[i];
      }
    }
  }
  return gradients;
}

Stencil Gradient::buildFaceNormalGradientStencil(const long& interId)
{

  bitpit::Interface inter = _grid->getInterface(interId);
  std::array<long,2> cellIds = inter.getOwnerNeigh();

  bool isUniform(true);
  int face = inter.getOwnerFace();
  int borderInterId;

  NPoint interCenter = _grid->evalInterfaceCentroid(interId);

  Stencil stencil;
  MLS2D interp;
  bool hasVertexOnBord(_grid->hasVertexOnBorder(interId, borderInterId));

  double h = _grid->evalInterfaceArea(interId);

  double sign(0.);
  switch(face%2)
  {
  case 0:
    sign = -1.;
    break;
  case 1:
    sign = 1.;
    break;
  default:
    break;
  }

  if (!hasVertexOnBord)
  {
    isUniform = ( ( _grid->getCellLevel(cellIds[0]) ==
                    _grid->getCellLevel(cellIds[1]) ) ? true : false);

    if (isUniform)
    {
      stencil.neighs.push_back(cellIds[0]);
      stencil.weights.weights.push_back(-1./h);
      stencil.isOutsideCell.push_back(0);

      stencil.neighs.push_back(cellIds[1]);
      stencil.weights.weights.push_back(1./h);
      stencil.isOutsideCell.push_back(0);

    }
    else
    {
      std::vector<long> neighs;
      std::vector<int>  faceVertex =
        _grid->getFaceVertexLocalIds(inter.getOwnerFace());

      _grid->findCellVertexOneRing(cellIds[0],
                                   faceVertex[0],
                                   &neighs);
      _grid->findCellVertexOneRing(cellIds[0],
                                   faceVertex[1],
                                   &neighs);

      std::vector<NPoint> xref;
      for (auto &nId: neighs)
      {
        xref.push_back(_grid->evalCellCentroid(nId));
      }

      std::vector<double> weights =
        interp.computeInterpolationAtFaceCenter(interCenter,
                                                xref);
      size_t nNeighs = neighs.size();

      // This offset tells us if we need to take
      // dx or dy derivative for normal gradient
      // Face 0 or 1: normal gradient = (-+) dx
      // Face 2 or 3: normal gradient = (-+) dy
      int offset = face/2 + 1;
      for (size_t i=0; i<nNeighs; i++)
      {
        stencil.neighs.push_back(neighs.at(i));
        stencil.weights.weights.push_back(sign * weights.at(i * 3 + offset));
        stencil.isOutsideCell.push_back(0);

      }
    }

  }
  else   // At least one vertex on boundary
  {
    if ( !(inter.isBorder()) )
    {
      if( _grid->getCellLevel(inter.getOwner()) ==
          _grid->getCellLevel(inter.getNeigh()) )
      {
        stencil.neighs.push_back(cellIds[0]);
        stencil.weights.weights.push_back(-1./h);
        stencil.isOutsideCell.push_back(0);

        stencil.neighs.push_back(cellIds[1]);
        stencil.weights.weights.push_back(1./h);
        stencil.isOutsideCell.push_back(0);
      }
      else
      {
        std::vector<long> neighs;
        std::vector<int>  faceVertex =
          _grid->getFaceVertexLocalIds(inter.getOwnerFace());

        _grid->findCellVertexOneRing(cellIds[0],
                                     faceVertex[0],
                                     &neighs);
        _grid->findCellVertexOneRing(cellIds[0],
                                     faceVertex[1],
                                     &neighs);

        std::vector<NPoint> xref;
        std::vector<int> isGhostCell;
        for (auto &nId: neighs)
        {
          xref.push_back(_grid->evalCellCentroid(nId));
          isGhostCell.push_back(0);
        }

        // Add now the two cells for boundary contribution
        NPoint xBuff(_grid->evalCellCentroid(inter.getOwner()));
        xBuff[borderInterId/2] += pow(-1, borderInterId+1)
                                  * _grid->evalCellSize(inter.getOwner());

        xref.push_back(xBuff);
        neighs.push_back(cellIds[0]);
        isGhostCell.push_back(-borderInterId - 1);

        xBuff = _grid->evalCellCentroid(inter.getNeigh());
        xBuff[borderInterId/2] += pow(-1, borderInterId+1)
                                  * _grid->evalCellSize(inter.getNeigh());

        xref.push_back(xBuff);
        neighs.push_back(cellIds[1]);
        isGhostCell.push_back(-borderInterId - 1);

        std::vector<double> weights =
          interp.computeInterpolationAtFaceCenter(interCenter,
                                                  xref);
        int nNeighs = neighs.size();
        int offset = face/2 + 1;
        for (int i=0; i<nNeighs; i++)
        {
          stencil.neighs.push_back(neighs.at(i));
          stencil.weights.weights.push_back(sign * weights.at(i * 3 + offset));
          stencil.isOutsideCell.push_back(isGhostCell.at(i));
        }
      }
    }
    else
    {
      double h = _grid->evalInterfaceArea(interId);
      stencil.neighs.push_back(cellIds[0]);
      stencil.weights.weights.push_back(-1./h);
      stencil.isOutsideCell.push_back(0);

      stencil.neighs.push_back(cellIds[0]);
      stencil.weights.weights.push_back(1./h);
      stencil.isOutsideCell.push_back(-borderInterId -1);
    }
  }
  return stencil;
}

PiercedVector<std::array<double, 3> > Gradient::computeFVGradient(
  PiercedVector<double>&  val,
  const double t,
  const Var type)
{

  // FIXME: HAS TO BE IMPLEMENTED FOR NON-UNIFORM CONFIGURATION
  PiercedVector<std::array<double,3> > gradient;
  Stencil stencil;

  if (Gradient::_userComm_only_one_init_dirty!=10101)
  {
    Gradient::_userComm = new UserDataComm<double>(*_grid, 1);
    Gradient::_userComm_only_one_init_dirty=10101;
  }

  Gradient::_userComm->update(1);
  Gradient::_userComm->communicate(val);


  auto& stencils = (_stencil->getInterfaceGradientStencil());
  Stencil sten;

  for (auto &inter: _grid->getInterfaces())
  {
    const long& interId = inter.getId();
    const auto& owners  = inter.getOwnerNeigh();
    const auto& normal = _grid->evalInterfaceNormal(interId);
    const double& area = _grid->evalInterfaceArea(interId);
    const double& volume = _grid->evalCellVolume(owners[0]);
    const int& face = inter.getOwnerFace();
    const int dim = face/2;

    if (_hasStencil)
      sten = stencils[interId];
    else
      sten = this->buildFaceNormalGradientStencil(interId);

    if (!inter.isBorder())
    {
      if (!gradient.exists(owners[0]))
        gradient.emplace(owners[0]);
      if (!gradient.exists(owners[1]))
        gradient.emplace(owners[1]);

      gradient[owners[0]][dim] +=
        0.5 * (val[owners[0]] + val[owners[1]] )
        * area/volume
        *std::accumulate(normal.begin(),
                         normal.end(),
                         0.);

      gradient[owners[1]][dim] -=
        0.5 * (val[owners[0]] + val[owners[1]] )
        * area/volume
        *std::accumulate(normal.begin(),
                         normal.end(),
                         0.);
    }
    else
    {
      if (!gradient.exists(owners[0]))
        gradient.emplace(owners[0]);
      NPoint ghostOctCenter =  _grid->evalCellCentroid(owners[0]);
      ghostOctCenter[dim] += pow(-1,face+1)
                             * _grid->evalCellSize(owners[0]);

      double value(0.);
      value = _grid->computeBoundaryValue(ghostOctCenter,
                                          val[owners[0]],
                                          face,
                                          type,
                                          t);
      gradient[owners[0]][dim] += value * area/volume
                                  *std::accumulate(normal.begin(),
                                                   normal.end(),
                                                   0.);
    }
  }

  //FIXME: Maybe useless ?
  _userComm->update(3);
  _userComm->communicate(gradient);
  return gradient;
}
}
