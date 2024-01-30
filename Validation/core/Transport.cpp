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
 * @file Transport.cpp
 * @brief This file contains the Transport class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 * @copyright Inria
 */


#include "Transport.hpp"
#include "Utils.hpp"
#include "InterpolatorFactory.hpp"
#include "Gradient.hpp"
#include "NeosAssert.hpp"
#include <numeric>
#include <chrono>

namespace neos {

UserDataComm<double>* Transport::_pvComm;
int Transport::_pvComm_only_one_init_dirty=0;

Transport::Transport(Grid* grid, StencilBuilder* stencil) : _grid(grid),
  _stencil(stencil)
{
}

double Transport::computeInterpo(PiercedVector<double> &data,
                                 const long &id,
                                 NPoint x_old)
{
  double res = 0;
  long id_old = _grid->locatePoint(x_old);
  _interpo = InterpolatorFactory::get(_interpType);
  if (_interpType == interpoType::RBF)
    ((Rbf*)_interpo)->setEpsilon(_epsilon);

  if (id_old >= 0) {
    std::vector<long> neighs = _grid->findQuadrantNeighs(x_old, id_old);
    std::vector<NPoint> x_ref(neighs.size());
    std::vector<double> phi_ref(neighs.size());


    for (std::size_t i = 0; i < neighs.size(); ++i) {
      x_ref[i] = _grid->evalCellCentroid(neighs[i]);
      phi_ref[i] = data[neighs[i]];
    }
    res = _interpo->computeInterpolation(x_old, x_ref, phi_ref);
  }
  else
  {
    res = data[id];
  }
  delete _interpo;
  return res;
}

void Transport::compute(PiercedVector<double> &data, const std::vector<double> &u, const double dt) {
  PiercedVector<double> ustar;
  ustar.reserve(_grid->nbCells());

  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior())
    {
      const long &id = cell.getId();
      if (!ustar.exists(id)) {
        ustar.emplace(id);
      }
      NPoint x0 = _grid->evalCellCentroid(id);
      NPoint x_old;
      for (int i = 0; i < _grid->getDim(); ++i) {
        x_old[i] = x0[i] - dt * u[i];
      }
      ustar[id] = computeInterpo(data, id, x_old);
    }
  }
  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior())
    {
      const long &id = cell.getId();
      data[id] = ustar[id];
    }
  }
}

void Transport::compute(PiercedVector<double> &data, const PiercedVector<NPoint> &u, const double dt) {
  PiercedVector<double> ustar;
  ustar.reserve(_grid->getCellCount());

  for (auto &cell : _grid->getCells()) {
    const long &id = cell.getId();

    if (!ustar.exists(id)) {
      ustar.emplace(id);
    }
    NPoint x0 = _grid->evalCellCentroid(id);
    NPoint x_old;

    x_old[NPX] = x0[NPX] - dt * u[id][NPX];
    x_old[NPY] = x0[NPY] - dt * u[id][NPY];
    x_old[NPZ] = x0[NPZ] - dt * u[id][NPZ];

    ustar[id] = computeInterpo(data, id, x_old);
  }
  for (auto &cell : _grid->getCells()) {
    const long &id = cell.getId();
    data[id] = ustar[id];
  }
}

std::vector<PiercedVector<double> > Transport::computeWithSecondOrderFV(
  const std::vector<PiercedVector<double> > &velFC,
  std::vector<PiercedVector<double> > &u,
  const std::vector<Var>& types,
  const std::vector<double>& times)
{


  // Assert than each vector has the same sizes
  // u: vector of transported variables
  // velFC: vector of centered velocities
  // types: vector of types for boundaries conditions
  // times: vector which contains current times for transported values
  // FIXME: Has to be encapsulate in a ifdef debug
  ASSERT(types.size() == u.size(),
         "Vector of types for boundary conditions has not the same size as \
    vector of transported values in Transport::computeWithSecondOrderFV \
    method. You have to specify a type of boundary conditions for each \
    transported value");
  ASSERT(velFC.size() == u.size(),
         "Vector of face centered velocities has not the same size as vector \
	   of transported values in Transport::computeWithSecondOrderFV \
	   method. You have to specify a velocity vector for each transported \
	   value");
  ASSERT((times.size() == 0) || (times.size() == u.size()),
         "Vector of times has not the same size as vector of transported \
	   values in Transport::computeWithSecondOrderFV method. You have to \
	   specify a time for each transported value");


  // Define all useful variables
  std::array<long, 2> owners;
  bool boundary;
  NPoint interCenter, octCenter, deltaX, normal;
  long interId;
  Gradient grad(_grid->getDimension(), _grid, _stencil);
  double uCenter0, uCenter1;

  size_t nbVal = u.size();
  int boundFace;
  std::vector<PiercedVector<std::array<double,3> > > gradients(nbVal);
  std::vector<PiercedVector<double> > fluxFC(nbVal);
  std::vector<PiercedVector<double> > divergences(nbVal);

  for (size_t i=0; i<nbVal; i++)
  {
    fluxFC[i].reserve(_grid->getInterfaceCount());
    divergences[i].reserve(_grid->nbCells());
  }

  // Compute Cell Centered Least Square gradients for
  // each transported values
  gradients = grad.computeCCLSGradient(u, types, times);

  if (Transport::_pvComm_only_one_init_dirty!=10101)
  {
    Transport::_pvComm = new UserDataComm<double>(*_grid, 3);
    Transport::_pvComm_only_one_init_dirty=10101;
  }

  // Communicate each gradient
  for (size_t i=0; i<nbVal; i++)
  {
    Transport::_pvComm->update();
    Transport::_pvComm->communicate(gradients[i]);
  }

  // loop over interfaces
  for (auto &inter: _grid->getInterfaces())
  {

    interId     = inter.getId();
    owners      = inter.getOwnerNeigh();
    boundary    = inter.isBorder();
    interCenter = _grid->evalInterfaceCentroid(interId);
    normal = _grid->evalInterfaceNormal(interId);

    for (size_t i=0; i<nbVal; i++)
    {
      double t = times[i];

      octCenter = _grid->evalCellCentroid(owners[0]);
      deltaX = interCenter  - octCenter;
      fluxFC[i].emplace(interId);

      // Reconstruct u at face center using Taylor expansion for both
      // owners
      uCenter0 = u[i][owners[0]]
                 + Gradient::DotProduct(gradients[i][owners[0]], deltaX);

      // Reconstruction of u for owners[1] depends if the interface
      // is on the boundary. If not, we can use u[owners[1]]. Else,
      // this values doesn't exist and we have to use boundary conditions.
      if(!boundary)
      {
        octCenter = _grid->evalCellCentroid(owners[1]);
        deltaX = interCenter  - octCenter;
        uCenter1 = u[i][owners[1]]
                   + Gradient::DotProduct(gradients[i][owners[1]], deltaX);
      }
      else
      {
        boundFace = inter.getOwnerFace();

        // u from owners[1] depends on the variable type.
        // Var::LS is for a level set function
        if (types[i] == Var::LS)
        {
          uCenter1 = uCenter0;
        }
        else
        {
          uCenter1 = _grid->computeBoundaryValue(interCenter,
                                                 u[i][owners[0]],
                                                 boundFace,
                                                 types[i],
                                                 t);
        }
      }

      // Compute the flux at face center using a Local Lax-Friedrich
      // scheme. See Antoine Fondaneche's Transport Survey.

      fluxFC[i][interId] = 0.5 * velFC[i][interId] * (uCenter0 + uCenter1)
                           - 0.5 * std::fabs(velFC[i][interId])
                           * (uCenter1 - uCenter0) * std::accumulate(normal.begin(),
                                                                     normal.end(),
                                                                     0);

      // Compute Divergence of the flux.
      // For each cell:
      // div(F_i) = 1/|\omega_i| \sum_{k=1}^{n_f} F_{k,fc} * |f_k|
      // F_i: flux for cell i
      // |\omega_i|: volume of cell i
      // |f_k|: area of interface k
      // F_{k,fc}: value of F at face k center.
      // n_f: number of face

      if (_grid->getCell(owners[0]).isInterior())
      {
        if (!divergences[i].exists(owners[0]))
          divergences[i].emplace(owners[0]);

        divergences[i][owners[0]] += fluxFC[i][interId]
                                     * std::accumulate(normal.begin(), normal.end(), 0)
                                     * _grid->evalInterfaceArea(interId)
                                     / _grid->evalCellVolume(owners[0]);
      }

      if (owners[1] != bitpit::Element::NULL_ID)
      {
        if (_grid->getCell(owners[1]).isInterior())
        {
          if (!divergences[i].exists(owners[1]))
            divergences[i].emplace(owners[1]);

          divergences[i][owners[1]] -= fluxFC[i][interId]
                                       * std::accumulate(normal.begin(), normal.end(), 0)
                                       * _grid->evalInterfaceArea(interId)
                                       / _grid->evalCellVolume(owners[1]);
        }
      }
    }
  }
  return divergences;
}

// FIXME: There is surely a best option to implement this (template ? )
PiercedVector<double> Transport::computeWithSecondOrderFV(
  const PiercedVector<double> &velFC,
  PiercedVector<double> &u,
  const double t,
  const Var type)
{
  // FIXME: Find a way to do that in another way.
  auto uVector = std::vector<PiercedVector<double> >{u};

  auto divergence = this->computeWithSecondOrderFV(
    std::vector<PiercedVector<double> >{velFC},
    uVector,
    std::vector<Var>{type},
    std::vector<double>{t});


  //PiercedVector<double> res;
  //for (auto &id: divergence[0].getIds())
  //{
  //	res.emplace(id);
  //	res[id] = divergence[0][id];
  //}

  return divergence[0];
}
}
