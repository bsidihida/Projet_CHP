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

from mpi4py import MPI
from pyneos import pyneos
import sys
import math
import numpy as np
from pprint import pprint

import inspect

level = 6
dimension= 2
time_iteration= 20

g = pyneos.Grid(0, 0, 0, 1, 1.0/pow(2,level), dimension)

g.setExportName("time-poisson");

lap = pyneos.LaplacianFactory.get(pyneos.lapType.FINITEVOLUME, pyneos.solverType.PETSC, g);

lap.setBCPositionCoef(pyneos.BCPosition.OnInterface);
lap.buildMatrix();
lap.setDirichletCondition();

internal =  g.internal()
afafb = [0] * g.nbCells()
vg0 = [0] * g.nbCells()
vg1 = [0] * g.nbCells()

#item = internal.getVTK()
#methodList = [method for method in dir(item) if callable(getattr(item, method))]
#processFunc = 1 and (lambda s: " ".join(s.split())) or (lambda s: s)
#print "\n".join(["%s %s" % (method.ljust(10), processFunc(str(getattr(item,
#                                                                       method).__doc__)))
#                  for method in methodList])

for i in xrange(time_iteration):
    for cell in internal.getCells():
        id =  cell.getId()
        if cell.isInterior():
            ic= internal.evalCellCentroid(id)
            x= ic[0]
            y= ic[1]
            g0= math.exp(-10*( pow(x+(i/(time_iteration*2.0))-0.6,2)+pow(y-(i/(time_iteration*2.0))-0.6,2))) + 0.001;
            g1= math.exp(-10*( pow(x-0.3,2)+pow(y-0.3,2) )) + 0.001;
            afafb[id]= (g0-g1)/g1
            vg0[id]= g0
            vg1[id]= g1

    lap.setRHS( afafb )
    lap.solveLaplacian()
    lap.addSolutionToVTK("solution")
    g.write()
