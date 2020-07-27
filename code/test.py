import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from DomainTriangulation import DomainTriangulation
from LiftTriangles import LiftTriangles
from PS_Spline import PS_Spline

# Prima triangolazione: stella a 4 punte.
points_x = np.array([8, 9.06, 12, 9.06, 8, 6.94, 4, 6.94, 8])
points_y = np.array([2, 4.94, 6, 7.06, 10, 7.06, 6, 4.94, 6])
triangles = np.array([[0, 1, 8],
                      [0, 7, 8],
                      [1, 2, 8],
                      [2, 3, 8],
                      [3, 4, 8],
                      [4, 5, 8],
                      [5, 6, 8],
                      [6, 7, 8]])
tr1 = Triangulation(points_x, points_y, triangles)

# Seconda triangolazione: stella a 5 punte.
points_x = np.array([4.95, 7.86, 10.79, 9.66, 12.57, 8.97, 7.84, 6.74, 3.13, 6.05, 7.86])
points_y = np.array([3.09, 5.22, 3.11, 6.53, 8.66, 8.65, 12.08, 8.64, 8.63, 6.52, 7.11])
tr = np.array([[0, 1, 9],
               [1, 2, 3],
               [3, 4, 5],
               [5, 6, 7],
               [7, 8, 9],
               [1, 9, 10],
               [1, 3, 10],
               [3, 5, 10],
               [5, 7, 10],
               [7, 9, 10]])
tr2 = Triangulation(points_x, points_y,tr)

# Terza triangolazione: esagono.
points_x = np.array([1, 2, 2.5, 2, 1, 0.5, 1.5])
points_y = np.array([0, 0, 0.87, 1.73, 1.73, 0.87, 0.87])
tr3 = Triangulation(points_x, points_y)

# Quarta triangolazione: freccia.
points_x = np.array([1, 3, 5, 5, 8, 5, 5, 1])
points_y = np.array([0, 0.5, 0, -1, 0.5, 2, 1, 1])
tr = np.array([[0, 2, 1],
               [1, 2, 6],
               [1, 6, 7],
               [0, 1, 7],
               [2, 3, 4],
               [2, 4, 6],
               [4, 5, 6]])
tr4 = Triangulation(points_x, points_y, tr)

# Quinta triangolazione: triangolo equilatero.
points_x = np.array([1, 2, 1.5])
points_y = np.array([0, 0, 0.87])
tr5 = Triangulation(points_x, points_y)

# Sesta triangolazione: due triangoli equilateri.
points_x = np.array([1, 2, 1.5, 0.5])
points_y = np.array([0, 0, 0.87, 0.87])
tr = np.array([[0,1,2],
               [0,2,3]])
tr6 = Triangulation(points_x, points_y,tr)

# ***********************************************************
# Figura: elemento di base relativo al vertice, 1.
'''
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
plt.text(lf.domain.startPoints[0,0]+incr,lf.domain.startPoints[0,1]-incr,'$b_{1,r}^v$',fontsize=8)
plt.plot(lf.domain.startPoints[0,0],lf.domain.startPoints[0,1],marker='o',color='r')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{9,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{10,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{11,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{12,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{13,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{14,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{15,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{16,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
'''

# ***********************************************************
# Figura: elemento di base relativo al vertice, 2.
"""
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
plt.text(lf.domain.startPoints[0,0]+incr,lf.domain.startPoints[0,1]-incr,'$b_{1,r}^v$',fontsize=8)
plt.plot(lf.domain.startPoints[0,0],lf.domain.startPoints[0,1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{9,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{10,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{11,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{12,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{13,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{14,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{15,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{16,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
lf.TSoll.mask = [False,True,True]
plt.triplot(lf.TSoll)
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
"""

# ***********************************************************
# Figura: elemento di base relativo al vertice, 3.
"""
lf = LiftTriangles(tr6)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
plt.text(lf.domain.startPoints[0,0]+incr,lf.domain.startPoints[0,1]-incr,'$b_{1,r}^v$',fontsize=8)
plt.plot(lf.domain.startPoints[0,0],lf.domain.startPoints[0,1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{9,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{10,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{11,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{12,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{13,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{14,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{15,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{16,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[1]
plt.text(p[0]+incr,p[1]-incr,'$\widetilde{b_{6,r}^v}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
plt.plot(lf.ESollVx[0:2,0],lf.ESollVx[0:2,1],color='#fac900')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
"""

# ***********************************************************
# Figura: elemento di base relativo al vertice, 4.
"""
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
plt.text(lf.domain.startPoints[0,0]+incr,lf.domain.startPoints[0,1]-incr,'$b_{1,r}^v$',fontsize=8)
plt.plot(lf.domain.startPoints[0,0],lf.domain.startPoints[0,1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{9,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{10,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 2/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{11,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{12,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{13,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{14,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{15,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{16,r}^v$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 2/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.startPoints[1]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.startPoints[1] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p =1/3*lf.domain.startPoints[1] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.startPoints[2]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.startPoints[2] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p =1/3*lf.domain.startPoints[2] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
"""

# ***********************************************************
# Figura: plotting elementi di base al vertice.
'''
varray = np.zeros((7,3))
varray[6,0]=1
tarray = np.zeros(36)
edict = {}
ps = PS_Spline(tr3, varray, tarray, edict)
ps.spline3DPlot(save=True)
ps.splineContourPlot(showLift=True, save=True)
varray[6,0]=0
varray[6,1]=1
ps = PS_Spline(tr3, varray, tarray, edict)
ps.spline3DPlot(save=True)
ps.splineContourPlot(showLift=True, save=True)
varray[6,1]=0
varray[6,2]=1
ps = PS_Spline(tr3, varray, tarray, edict)
ps.spline3DPlot(save=True)
ps.splineContourPlot(showLift=True, save=True)
'''


# ***********************************************************
# Figura: elemento di base relativo al triangolo PS, 1.
'''
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{1}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[2] + 1/3*lf.domain.startPoints[1] + 1/3*lf.domain.startPoints[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.TSollVx[lf.TSollTr[0,2]]
plt.text(p[0]+incr,p[1],'$1$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='black')
p = lf.TSollVx[lf.TSollTr[0,1]]
plt.text(p[0]+incr,p[1],'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='black')
p = lf.TSollVx[lf.TSollTr[0,0]]
plt.text(p[0]+incr,p[1],'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='black')
lf.TSoll.mask = [False,True,True]
plt.triplot(lf.TSoll)
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
'''

# ***********************************************************
# Figura: elemento di base relativo al triangolo PS, 2.
'''
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{1}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[2] + 1/3*lf.domain.startPoints[1] + 1/3*lf.domain.startPoints[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[1] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[2] + 1/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[2] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
'''

# ***********************************************************
# Figura: elemento di base relativo al triangolo PS, 3.
'''
lf = LiftTriangles(tr6)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{1}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[2] + 1/3*lf.domain.startPoints[1] + 1/3*lf.domain.startPoints[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8}^t$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='b')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.incenters[1]+1/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
p = 1/3*lf.domain.incenters[1]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$0$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='g')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
'''

# ***********************************************************
# Figura: elemento di base relativo al lato di bordo.
'''
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{1}^e$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
p = lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2}^e$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='r')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
'''

# ***********************************************************
# Figura: plotting elementi di base al triangolo e bordo.
'''
varray = np.zeros((7,3))
tarray = np.zeros(36)
tarray[0] = 1
edict = {}
ps = PS_Spline(tr3, varray, tarray, edict)
ps.spline3DPlot(save=True)
ps.splineContourPlot(showLift=True, save=True)
varray = np.zeros((7,3))
tarray = np.zeros(36)
edict = {(0,1):1}
ps = PS_Spline(tr3, varray, tarray, edict)
ps.spline3DPlot(save=True)
ps.splineContourPlot(showLift=True, save=True)
'''

# ***********************************************************
# Figura: plotting pandoro.
'''
varray = np.zeros((7,3))
varray[6] = 1
tarray = np.ones(36)
edict = {}
ps = PS_Spline(tr3, varray, tarray, edict)
ps.spline3DPlot(save=True,cmap='RdBu')
ps.splineContourPlot(save=True, cmap='RdBu')
'''


# ***********************************************************
# Figura: plotting stella a 5 punte.
'''
varray = np.zeros((11,3))
tarray = np.ones(60)
edict = {}
ps = PS_Spline(tr2, varray, tarray, edict)
ps.spline3DPlot(save=True,cmap='Spectral',axIncl=(60,30))
ps.splineContourPlot(save=True, cmap='Spectral')
'''


# ***********************************************************
# Figura: plotting stella a 4 punte.
'''
varray = np.zeros((9,3))
tarray = np.zeros(48)
tarray[44] = 1
tarray[8] = 1
varray[7] = 1
tarray[40] = 1
tarray[32] = 1
varray[5] = 1
tarray[20] = 1
tarray[28] = 1
varray[3] = 1
tarray[2] = 1
tarray[16] = 1
varray[1] = 1
varray[8] = -2
edict = {(0,1):1,(1,0):1,(1,2):1,(2,1):1,(3,2):1,(2,3):1,(3,4):1,(4,3):1,
          (4,5):1,(5,4):1,(5,6):1,(6,5):1,(6,7):1,(7,6):1,(7,0):1,(0,7):1}
ps = PS_Spline(tr1, varray, tarray, edict)
ps.spline3DPlot(save=True,lvs=120)
ps.splineContourPlot(save=True)
'''

# ***********************************************************
# Figura: plotting freccia.
'''
varray = np.zeros((8,3))
varray[1] = 1
tarray = np.zeros(42)
tarray[24:] = 1
tarray[2] = 1
tarray[5] = 1
tarray[12] = 1
tarray[16] = 1
tarray[6] = 1
tarray[10] = 1
tarray[19] = 1
tarray[20] = 1
tarray[7:10] = 1
tarray[11] = 1
earray={}
ps = PS_Spline(tr4, varray, tarray, earray)
ps.spline3DPlot(save=True, cmap='winter', axIncl=(60,30))
ps.splineContourPlot(save=True, cmap='winter')
'''


# ***********************************************************
# Figura: insieme di determinazione esagono.
'''
dm = DomainTriangulation(tr3)
dm.domainPlot(save=True, showDP=True)
'''

# ***********************************************************
# Figura: altezze di bezier cumulate.
'''
lf = LiftTriangles(tr5)
lf.domain.domainPlot(showDP=True, colorPS='#5a5e99', legend=False, showExecute=False)
incr= 0.023
plt.axis('equal')
plt.text(lf.domain.startPoints[0,0]+incr,lf.domain.startPoints[0,1]-incr,'$b_{1}$',fontsize=8)
plt.plot(lf.domain.startPoints[0,0],lf.domain.startPoints[0,1],marker='o',color='grey')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{2}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{3}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 2/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{4}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{5}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,1][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{6}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{7}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.startPoints[0]+1/3*lf.domain.splitPoints[0,2][0]+1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{8}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.startPoints[0]+2/3*lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{9}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = lf.domain.splitPoints[0,1][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{10}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 2/3*lf.domain.splitPoints[0,1][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{11}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.splitPoints[0,1][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{12}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{13}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 1/3*lf.domain.splitPoints[0,2][0] + 2/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{14}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = 2/3*lf.domain.splitPoints[0,2][0] + 1/3*lf.domain.incenters[0]
plt.text(p[0]+incr,p[1]-incr,'$b_{15}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
p = lf.domain.splitPoints[0,2][0]
plt.text(p[0]+incr,p[1]-incr,'$b_{16}$',fontsize=8)
plt.plot(p[0],p[1],marker='o',color='grey')
nome = input('Nome del file?\n')
plt.savefig(nome, dpi=600)
plt.show()
'''

varray = np.zeros((11,3))
varray[10] = 1
tarray = np.ones(60)
edict = {}
ps = PS_Spline(tr2, varray, tarray, edict)
ps.splineContourPlot(save=True)












