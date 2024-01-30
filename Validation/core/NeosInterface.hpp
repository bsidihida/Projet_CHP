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

#ifndef __NEOSINT_HPP__
#define __NEOSINT_HPP__

/**
 * @file NeosInterface.hpp
 * @brief Neos insterface
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <Neos.hpp>

namespace neos {

class NeosInterface {
private:
Grid *_grid = NULL;
LevelSet *_ls = NULL;
Laplacian *_ilp = NULL;
uint8_t _dim = GRID_3D;

public:
NeosInterface() {
  bitpit::log::manager().initialize(bitpit::log::COMBINED);
  bitpit::log::cout() << consoleVerbosity(bitpit::log::QUIET);
}

~NeosInterface() {
  delete _grid;
  delete _ls;
}

void createGrid(NPoint origin, double L, double dh, uint8_t dim = GRID_3D);
int addSphere(NPoint center, double radius);
int addStl(const std::string &file, std::string tag);
int addSphere(NPoint center, double radius, std::string tag);
void write(std::string tag);
void write();
void addData(std::vector<double> data, std::string tag);
void addData(PiercedVector<double> &data, std::string tag);
void writeData();
void transport(std::string tag, const std::vector<double> &u, const double dt);
void transport(std::string tag, const std::vector<NPoint> &u, const double dt);
PiercedVector<double> laplacian(lapType lap, solverType solver, const PiercedVector<double> &rhs);
PiercedVector<double> laplacian(lapType lap, solverType solver, const std::vector<double> &rhs);
std::vector<NPoint> randomNPointVector();
std::vector<double> getLevelSet();
std::vector<double> getLevelSet(int id);
std::vector<double> getLevelSet(std::string tag);
void refine();
};
}

#endif /* __NEOSINT_HPP__ */
