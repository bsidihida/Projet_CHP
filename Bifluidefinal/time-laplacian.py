# -------------------------------------------------------------------------
#
#  NEOS
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of Neos.
#
#  Neos is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  Neos is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with Neos. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------

from pyneos.pyneos import neos
from mpi4py import MPI
import sys
import numpy as np
from pprint import pprint

#help(neos)
print(dir(neos.Grid))

level = 4
dimension= 3
temperature_core = 0.8
temperature_border = 0.2
temperature_ambiant = 0.2
coef_penalization = 1.0e20
tolerance_penalization = 0.1
coef_refinment = 0.2
delta_t= 0.001
time_iteration= 90
time_enable= 1
border_enable= 0
h=1

if (sys.argv[1] == "1"):
    grid = neos.Grid(0, 0, 0, 1, 1.0/pow(2,level), dimension)
    geo = neos.ASphere(0.5, 0.5, 0.5, 0.2, dimension)
    geo.computeLevelSet(grid);
    ls = neos.LevelSet(grid);
    geoId = ls.addGeometry(geo)
else:
    grid = neos.Grid(-20, -20, -10, 50.0, 2, dimension)
    geo = neos.STLGeometry()
    geo.addSTL("./Dragonite.stl", grid)
    geo.computeLevelSet(grid)
    ls = neos.LevelSet(grid)
    geoId = ls.addGeometry(geo)
    delta_t= 0.1
    coef_refinment = 10.0

phi0 = ls.getLevelSet()

nbRef = grid.refineFromRangeVector(phi0,0.0,coef_refinment);

print(nbRef," Cells Refined")

geo.computeLevelSet(grid)
phi0 = ls.getLevelSet()
pen = np.array(phi0)
initial_temperature_RHS = np.array(phi0)

h_h= pow(h,dimension)
h_h_sur_delta_t= h_h/delta_t

for i in range(len(pen)):
  if (pen[i] > tolerance_penalization):
    pen[i] = 0.0
  else:
    pen[i] = 1.0

fixed_temperature_RHS = pen * temperature_core

for i in range(len(initial_temperature_RHS)):
    if initial_temperature_RHS[i] > tolerance_penalization:
        initial_temperature_RHS[i] = temperature_ambiant
    else:
        initial_temperature_RHS[i] = temperature_core

sol_v = np.array(initial_temperature_RHS)

grid.setExportName("time-laplacian")
grid.setCounter()

lap = neos.LaplacianFactory.get(neos.lapType.FINITEVOLUME, neos.solverType.PETSC, grid);

#lap.setBCPositionCoef(neos.BCPosition.OnInterface);
lap.buildMatrix();
lap.penalize(pen * coef_penalization * (-1.0) * h_h - time_enable * h_h_sur_delta_t);

for i in range(time_iteration):
    lap.setRHS( fixed_temperature_RHS * (coef_penalization * (-1.0) * h_h) -
               (time_enable * h_h_sur_delta_t) * sol_v)
    norm = lap.solveLaplacian()
    print("Iter ", i,  " Norm ", norm)
    sol_v = lap.getSolution()
    lap.addSolutionToVTK("solution");
    grid.addData("phi",phi0)
    grid.write()
    sol_v = np.array(sol_v)
