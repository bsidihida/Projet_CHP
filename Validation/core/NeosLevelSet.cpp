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
 * @file LevelSet.cpp
 * @brief This file contains the LevelSet meta classe
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <map>
#include "NeosLevelSet.hpp"
#include "Utils.hpp"

namespace neos {

int LevelSet::addGeometry(Geometry *geo, std::string tag) {
  std::map<std::string, int>::iterator itr;
  int uid;

  itr = _lvlstMap.find(tag);
  if (itr == _lvlstMap.end()) {
    ILevelSet *lvlst = new ILevelSet(_grid);
    uid = Utils::getUid();
    nbGeo++;
    lvlst->addGeometry(geo);
    _lvlst.insert(std::pair<int, ILevelSet*>(uid, lvlst));
    _lvlstMap.insert(std::pair<std::string, int>(tag, uid));
  } else {
    std::map<int, ILevelSet*>::iterator itGeo;
    uid = itr->second;
    itGeo = _lvlst.find(uid);
    itGeo->second->addGeometry(geo);
  }
  return uid;
}

int LevelSet::addGeometry(Geometry *geo) {
  ILevelSet *lvlst = new ILevelSet(_grid);
  int uid = Utils::getUid();

  nbGeo++;
  lvlst->addGeometry(geo);
  _lvlst.insert(std::pair<int, ILevelSet*>(uid, lvlst));
  return uid;
}

void LevelSet::delGeometry(int id) {
  std::map<std::string, int>::iterator itr = _lvlstMap.begin();

  nbGeo--;
  _lvlst.erase(id);
  while (itr != _lvlstMap.end()) {
    if (itr->second == id) {
      _lvlstMap.erase(itr++);
    } else {
      ++itr;
    }
  }
}

Geometry *LevelSet::getGeometry(int idGeo) {
  return _lvlst.find(idGeo)->second->getGeometry();
}

int LevelSet::getId(std::string tag) {
  return _lvlstMap.find(tag)->second;
}

void LevelSet::delGeometry(std::string tag) {
  std::map<std::string, int>::iterator itr;

  itr = _lvlstMap.find(tag);
  if (itr != _lvlstMap.end()) {
    nbGeo -= _lvlst.erase(itr->second);
    _lvlstMap.erase(itr);
  }
}

void LevelSet::printMaps() {
  for (std::map<std::string, int>::iterator itr = _lvlstMap.begin(); itr != _lvlstMap.end(); ++itr) {
    std::cout << itr->first << " " << itr->second << std::endl;
  }
  for (std::map<int, ILevelSet*>::iterator itr = _lvlst.begin(); itr != _lvlst.end(); ++itr) {
    std::cout << itr->first << " " << itr->second->getGeoNbr() << std::endl;
  }

}

double LevelSet::getLevelSet(const NPoint &cell) {
  double phi =  std::numeric_limits<double>::max();

  for (std::map<int, ILevelSet*>::iterator itr = _lvlst.begin(); itr != _lvlst.end(); ++itr) {
    double tphi = itr->second->getLevelSet(cell);
    if (tphi < phi) {
      phi = tphi;
    }
  }
  return phi;
}

double LevelSet::getLevelSet(const NPoint &cell, int id) {
  return _lvlst.find(id)->second->getLevelSet(cell);
}

double LevelSet::getLevelSet(const NPoint &cell, std::string tag) {
  int id = _lvlstMap.find(tag)->second;

  return _lvlst.find(id)->second->getLevelSet(cell);
}


std::vector<double> LevelSet::getLevelSet() {
  PiercedVector<double> fphi;

  fphi.reserve(_grid->nbCells());
  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior()) {
      const long &i=cell.getId();
      fphi.emplace(i);
    }
  }

  for (auto &cell : _grid->getCells()) {
    if (cell.isInterior()) {
      const long &id = cell.getId();
      double phi = std::numeric_limits<double>::max();
      for (std::map<int, ILevelSet*>::iterator itr = _lvlst.begin(); itr != _lvlst.end(); ++itr) {
        double tphi = itr->second->getLevelSet(id);
        if (tphi < phi) {
          phi = tphi;
        }
      }
      fphi[cell.getId()] = phi;
    }
  }

  std::vector<double> phi = PVtoV(fphi, _grid);
  return phi;
}

std::vector<double> LevelSet::getLevelSet(int id) {
  return _lvlst.find(id)->second->getLevelSet();
}

std::vector<double> LevelSet::getLevelSet(std::string tag) {
  int id = _lvlstMap.find(tag)->second;

  return _lvlst.find(id)->second->getLevelSet();
}

void LevelSet::transport(std::string tag, const PiercedVector<NPoint> &u, const double dt) {
  int id = _lvlstMap.find(tag)->second;

  _lvlst.find(id)->second->transport(_grid, u, dt);
}
}
