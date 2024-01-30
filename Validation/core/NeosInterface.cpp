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
 * @file NeosInterface.cpp
 * @brief Neos insterface
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <NeosInterface.hpp>
#include "LaplacianFactory.hpp"
#include "Grid.hpp"

namespace neos {

  void NeosInterface::createGrid(NPoint origin, double L, double dh, uint8_t dim) {
  _grid = new Grid(origin[NPX], origin[NPY], origin[NPZ], L, dh, dim);
  _grid->getVTK().setCounter(0);
  _dim = dim;
  _ls = new LevelSet(_grid);
}

int NeosInterface::addSphere(NPoint center, double radius) {
  ASphere       *geo = new ASphere(center[NPX], center[NPY], center[NPZ], radius, _dim);
  geo->computeLevelSet(_grid);
  return _ls->addGeometry(geo);
}

int NeosInterface::addStl(const std::string &file, std::string tag) {
  STLGeometry *geo = new STLGeometry(file, _grid);
  return _ls->addGeometry(geo, tag);
}

int NeosInterface::addSphere(NPoint center, double radius, std::string tag) {
  ASphere       *geo = new ASphere(center[NPX], center[NPY], center[NPZ], radius, _dim);
  geo->computeLevelSet(_grid);
  return _ls->addGeometry(geo, tag);
}

void NeosInterface::write(std::string tag) {
  std::vector<double> phi = _ls->getLevelSet(tag);

  _grid->getVTK().addData("phi", bitpit::VTKFieldType::SCALAR,
                          bitpit::VTKLocation::CELL, phi);
  _grid->write();
}

void NeosInterface::addData(PiercedVector<double> &data, std::string tag) {
  std::vector<double> vData = PVtoV(data, _grid);

  _grid->getVTK().addData(tag, bitpit::VTKFieldType::SCALAR,
                          bitpit::VTKLocation::CELL, vData);
}

void NeosInterface::addData(std::vector<double> data, std::string tag) {
  _grid->getVTK().addData(tag, bitpit::VTKFieldType::SCALAR,
                          bitpit::VTKLocation::CELL, data);
}

void NeosInterface::writeData() {
  _grid->write();
}

void NeosInterface::write() {
  std::vector<double> phi = _ls->getLevelSet();
  _grid->getVTK().addData("phi", bitpit::VTKFieldType::SCALAR,
                          bitpit::VTKLocation::CELL, phi);
  _grid->write();
  _grid->getVTK().removeData("phi");
}

void NeosInterface::transport(std::string tag, const std::vector<double> &u, const double dt) {
  PiercedVector<NPoint> PVu;

  for (auto &cell : _grid->getCells()) {
    const long &i = cell.getId();
    PVu.emplace(i);
    PVu[i] = { { u[NPX], u[NPY], u[NPZ] } };
  }
  _ls->transport(tag, PVu, dt);
}

void NeosInterface::transport(std::string tag, const std::vector<NPoint> &u, const double dt) {
  PiercedVector<NPoint> PVu;

  PVu = VtoPV(u, _grid);
  _ls->transport(tag, PVu, dt);
}

std::vector<NPoint> NeosInterface::randomNPointVector() {
  std::vector<NPoint> V;
  V.resize(_grid->getCellCount());

  for (size_t i = 0; i < V.size(); i++) {
    V[i][0] = 0.3;//sin(i);//(rand()/(double)RAND_MAX) * 1.0;
    V[i][1] = cos(i);//(rand()/(double)RAND_MAX) * 1.0;
    V[i][2] = 0.3;//cos(i);//tan(i);//(rand()/(double)RAND_MAX) * 1.0;
  }
  return V;
}

PiercedVector<double> NeosInterface::laplacian(lapType lap, solverType solver, const PiercedVector<double> &rhs) {
  _ilp = LaplacianFactory::get(lap, solver, _grid);
  _ilp->setBCPositionCoef(BCPosition::OnInterface);
  _ilp->setRHS(rhs);
  return _ilp->solve();
}

PiercedVector<double> NeosInterface::laplacian(lapType lap, solverType solver, const std::vector<double> &rhs) {
  _ilp = LaplacianFactory::get(lap, solver, _grid);
  _ilp->setBCPositionCoef(BCPosition::OnInterface);
  _ilp->setRHS(rhs);
  return _ilp->solve();
}

std::vector<double> NeosInterface::getLevelSet() {
  return _ls->getLevelSet();
}

std::vector<double> NeosInterface::getLevelSet(int id) {
  return _ls->getLevelSet(id);
}

std::vector<double> NeosInterface::getLevelSet(std::string tag) {
  return _ls->getLevelSet(tag);
}

void NeosInterface::refine() {
  for (auto &cell : _grid->getCells()) {
    const long &id = cell.getId();
    _grid->markCellForRefinement(id);
  }
  _grid->update();
}
}
