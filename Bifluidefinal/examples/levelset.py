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
import numpy as np

import inspect

#help(pyneos)
print(dir(neos))

dim = 3
# grid = pyneos.Grid(-20, -20, -20, 50.0, 0.5, dim) => 2 millions d'octans => crash
grid = neos.Grid(-20, -20, -20, 20, 0.5, dim)
ls= neos.LevelSet(grid);
#geo = pyneos.ASphere(10, 15, 10, 2.0, dim)
#geo2 = pyneos.ASphere(5, 7, 13, 2.0, dim)
#geo3 = pyneos.ASphere(7, 7, 3, 2.0, dim)
geo4 = neos.STLGeometry()
geo4.addSTL("./Dragonite.stl", grid);
ls= neos.LevelSet(grid);

geo.computeLevelSet(grid);
geo2.computeLevelSet(grid);
geo3.computeLevelSet(grid);

dir(grid)

ls.addGeometry(geo)

ls.addGeometry(geo, "TAG");
ls.addGeometry(geo2)
#ls.addGeometry(geo2, "TAG");
ls.addGeometry(geo4, "TAG");
ls.addGeometry(geo3,"TAGO")

print ("Number of geometry loaded : " , ls.countGeometry())
#ls.printMaps()

a = [0.0,0.0,0.0]
print (a)

i= ls.getLevelSet
print (i)
i= ls.getLevelSet(a,"TAG")
print (i)
i= ls.getLevelSet(a)
print (i)

# Cette methode n'existe pas pour la geometrie en question
j= ls.getLevelSet(2)
g= ls.getLevelSet("TAGO")


#g.printMaps()
