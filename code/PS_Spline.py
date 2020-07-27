import math as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.tri import Triangulation, UniformTriRefiner
from numpy.linalg import norm, solve
from mpl_toolkits.mplot3d import Axes3D
from DomainTriangulation import DomainTriangulation
from LiftTriangles import LiftTriangles

class PS_Spline:
    def __init__(self, tr: Triangulation, Vcoeffs: np.array, Tcoeffs: np.array, Ecoeffs: dict):
        # Inizializzazione oggetti legati al dominio.
        self.domain = DomainTriangulation(tr)
        self.lift_triangles = LiftTriangles(tr)

        # Controllo di conformita' di array e dizionario dei coefficienti.
        if Vcoeffs.shape == (len(self.domain.startPoints),3):
            self.Vcoeffs = Vcoeffs
        else:
            raise ValueError('Array dei coefficienti ai vertici non conforme: richiesto (' + \
                str(len(self.domain.startPoints)) + ',3).')

        if Tcoeffs.shape == (len(self.domain.PS_triangles),):
            self.Tcoeffs = Tcoeffs
        else:
            raise ValueError('Array dei coefficienti ai triangoli PS non conforme: richiesto ('+ \
                str(len(self.domain.PS_triangles)) + ',).')

        # Completamento del dizionario dei coefficienti ai lati di bordo.
        BoundEdgesList = []
        self.Ecoeffs = Ecoeffs
        for edge in self.domain.splitPoints:
            if self.domain.splitPoints[edge][1]:
                BoundEdgesList.append((edge[0],edge[1]))
                BoundEdgesList.append((edge[1],edge[0]))
        
        for el in BoundEdgesList:
            if not el in self.Ecoeffs:
                self.Ecoeffs[el] = 0

        # Calcolo delle altezze dei punti di controllo.
        self.setControlPoints()

        # Calcolo dei valori immagine della superficie.
        self.setSurface()

    def setControlPoints(self):
        # In controlPointsHeights[m,(i,j,k)] viene memorizzata l'altezza relativa al punto
        # di controllo sul triangolo PS m individuato dalla tripla (i,j,k). 
        self.controlPointsHeights = {}

        # Altezze sui vertici (punti (3,0,0))
        for i in range(len(self.domain.startPoints)):
            for tr,_,_ in self.domain.TList[i]:
                alpha = self.lift_triangles.alpha[i][0:3]
                self.controlPointsHeights[tr,(3,0,0)] = alpha[0]*self.Vcoeffs[i,0] + alpha[1]*self.Vcoeffs[i,1] + alpha[2]*self.Vcoeffs[i,2]

        # Primo round di calcolo delle altezze.
        for n in range(len(self.domain.startTr.triangles)):
            inHeight = []
            for i in self.domain.startTr.triangles[n]:
                # Ricavo la lista PSList dei triangoli PS nel triangolo n afferenti al
                # vertice i e dizionario jDict del corrispettivo j adiacente ad i.
                PSList, jList = self.lift_triangles.getPSTriangles(n,i)
                inHeight.append((self.domain.xi[n,i],PSList[0]))

                # Recupero i coefficienti alpha.
                alpha = self.lift_triangles.alpha[i]

                for m_index in range(len(PSList)):
                    m = PSList[m_index]
                    m_prime = PSList[(m_index+1)%2]
                    j = jList[m_index]
                    k = jList[(m_index+1)%2]

                    # punti (2,0,1)
                    (i_ord, j_ord) = (min(i,j),max(i,j))
                    diffr = self.domain.splitPoints[i_ord,j_ord][0] - self.domain.startPoints[i]
                    b = np.zeros(3)
                    b[0] = alpha[0] + 1/3*alpha[3]*diffr[0] + 1/3*alpha[6]*diffr[1]
                    b[1] = alpha[1] + 1/3*alpha[4]*diffr[0] + 1/3*alpha[7]*diffr[1]
                    b[2] = alpha[2] + 1/3*alpha[5]*diffr[0] + 1/3*alpha[8]*diffr[1]
                    self.controlPointsHeights[m,(2,0,1)] = b.dot(self.Vcoeffs[i])

                    # punti (2,1,0)
                    diffz = self.domain.incenters[n] - self.domain.startPoints[i]
                    b = np.zeros(3)
                    b[0] = alpha[0] + 1/3*alpha[3]*diffz[0] + 1/3*alpha[6]*diffz[1]
                    b[1] = alpha[1] + 1/3*alpha[4]*diffz[0] + 1/3*alpha[7]*diffz[1]
                    b[2] = alpha[2] + 1/3*alpha[5]*diffz[0] + 1/3*alpha[8]*diffz[1]
                    self.controlPointsHeights[m,(2,1,0)] = b.dot(self.Vcoeffs[i])

                    # punti (1,1,1)
                    b2 = self.domain.lamb[j,i]/(2*self.lift_triangles.TSollSigma[m])
                    b1 = 1-b2
                    self.controlPointsHeights[m,(1,1,1)] = b1*self.controlPointsHeights[m,(2,1,0)] + b2*self.Tcoeffs[m]
                    
                    # punti (1,2,0)
                    b2 = self.domain.xi[n,j]/(2*self.lift_triangles.TSollSigma[m])
                    b3 = self.domain.xi[n,k]/(2*self.lift_triangles.TSollSigma[m_prime])
                    b1 = 1 - b2 - b3
                    r = b1*self.controlPointsHeights[m,(2,1,0)] + b2*self.Tcoeffs[m] + b3*self.Tcoeffs[m_prime]
                    self.controlPointsHeights[m,(1,2,0)] = r

            # punti (0,3,0): altezze relative agli incentri.
            for i in self.domain.startTr.triangles[n]:
                PSList,_ = self.lift_triangles.getPSTriangles(n,i)
                for m in PSList:
                    b1 = self.controlPointsHeights[inHeight[0][1],(1,2,0)]
                    b2 = self.controlPointsHeights[inHeight[1][1],(1,2,0)]
                    b3 = self.controlPointsHeights[inHeight[2][1],(1,2,0)]
                    r = inHeight[0][0]*b1 + inHeight[1][0]*b2 + inHeight[2][0]*b3
                    self.controlPointsHeights[m,(0,3,0)] = r

        # Secondo round di calcolo delle altezze.
        for n in range(len(self.domain.startTr.triangles)):
            for i in self.domain.startTr.triangles[n]:
                # Ricavo la lista PSList dei triangoli PS nel triangolo n afferenti al
                # vertice i e dizionario jDict del corrispettivo j adiacente ad i.
                PSList, jList = self.lift_triangles.getPSTriangles(n,i)

                for m_index in range(len(PSList)):
                    m = PSList[m_index]
                    j = jList[m_index]

                    # controllo se il lato ViRij sta sul bordo del dominio
                    (i_ord,j_ord) = (min(i,j),max(i,j))
                    if self.domain.splitPoints[i_ord,j_ord][1]:
                        # punti (1,0,2), caso al bordo
                        b2 = self.domain.lamb[j,i]/(2*self.lift_triangles.ESollSigma[i,j])
                        b1 = 1-b2
                        r = b1*self.controlPointsHeights[m,(2,0,1)] + b2*self.Ecoeffs[i,j]
                        self.controlPointsHeights[m,(1,0,2)] = r
                    else:
                        # punti (1,0,2), caso non al bordo
                        n_prime = self.domain.getCommonTriangle(n,(i,j))
                        m_prime = self.domain.getRijOpposite(n,(i,j))
                        b1 = self.domain.mu[n,n_prime]
                        b2 = self.domain.mu[n_prime,n]
                        r = b1*self.controlPointsHeights[m,(1,1,1)] + b2*self.controlPointsHeights[m_prime,(1,1,1)]
                        self.controlPointsHeights[m,(1,0,2)] = r

        # Terzo ed ultimo round di calcolo delle altezze.
        for n in range(len(self.domain.startTr.triangles)):
            (PSList, ijList) = self.domain.getZnOpposites(n)

            for m_index in range(len(PSList)):
                (m,m_prime) = PSList[m_index]
                (i,j) = ijList[m_index]
                b1 = self.domain.lamb[i,j]
                b2 = self.domain.lamb[j,i]

                # punti (0,0,3)
                r = b1*self.controlPointsHeights[m,(1,0,2)] + b2*self.controlPointsHeights[m_prime,(1,0,2)]
                self.controlPointsHeights[m,(0,0,3)] = r

                # punti (0,1,2)
                r = b1*self.controlPointsHeights[m,(1,1,1)] + b2*self.controlPointsHeights[m_prime,(1,1,1)]
                self.controlPointsHeights[m,(0,1,2)] = r

                # punti (0,2,1)
                r = b1*self.controlPointsHeights[m,(1,2,0)] + b2*self.controlPointsHeights[m_prime,(1,2,0)]
                self.controlPointsHeights[m,(0,2,1)] = r

    def setSurface(self):
        # Crea i valori immagine per la superficie, calcolando su ogni punto della plotTriangulation il valore della superficie di Bezier triangolare sul triangolo PS al quale il punto appartiene.

        # Raffinamento uniforme della triangolazione Powell_Sabin ai fini del plottaggio.
        refiner = UniformTriRefiner(self.domain.PowellSabinTr)
        self.plotTriangulation, self.TrIndex = refiner.refine_triangulation(return_tri_index=True, subdiv=4)

        # Salvo i punti della triangolazione iniziale in un solo array.
        self.plotPoints = np.array(list(zip(self.plotTriangulation.x,self.plotTriangulation.y)))

        # Inizializzo l'array delle ordinate.
        plotHeights = []

        for point_index in range(len(self.plotPoints)):
            point = self.plotPoints[point_index]

            # Ricerco il triangolo PS al quale appartiene il punto.
            m = self.TrIndex[point_index]
            triangle = self.domain.PS_triangles[m]
            tr_vertexes = np.array([self.domain.PS_points[i] for i in triangle])

            # Ricavo le coordinate baricentriche bc del punto rispetto al triangolo PS.
            bc = self.getBarycentricCoords(tr_vertexes,point)

            # Calcolo l'ordinata associata al punto.
            plotHeights.append(self.triDeCasteljau(m,(0,0,0),bc,3))

            # Salvo array delle ordinate.
            self.plotHeights = np.array(plotHeights)

    def getBarycentricCoords(self, triangle: np.array, point: np.array):
        # Ricava le coordinate baricentriche rispetto a triangle del punto point.
        A = np.block([[triangle[:,0]],[triangle[:,1]],[1,1,1]])
        b = np.array([point[0], point[1], 1])
        return solve(A,b)

    def triDeCasteljau(self, m:int, index:tuple, bc:np.array, deg:int):
        # Esegue ricosivamente l'algoritmo di DeCasteljau triangolare per il calcolo
        # dell'ordinata relativa al punto di coordinate baricentriche bc rispetto
        # al triangolo m.
        if deg == 0:
            return self.controlPointsHeights[m,index]
        else:
            i0 = (index[0]+1,index[1],index[2])
            i1 = (index[0],index[1]+1,index[2])
            i2 = (index[0],index[1],index[2]+1)
            return bc[0]*self.triDeCasteljau(m,i0,bc,deg-1) + bc[1]*self.triDeCasteljau(m,i1,bc,deg-1) + bc[2]*self.triDeCasteljau(m,i2,bc,deg-1)


    def spline3DPlot(self, cmap='twilight_shifted', curveLevels=True, axis=False, save=False, lvs=60, axIncl=(40,270)):
        # Visualizza la superficie B-spline.
        # cmap: la mappa colori da utilizzare sulle linee di livello;
        # cureLevels: se True visualizza la superficie sottoforma di curve-livello;
        # axis: se True mostra gli assi cartesiani;
        # save: se True salva l'immagine sottoforma di file;
        # lvs: numero di livelli da mostrare;
        # axIncl: posizione iniziale degli assi. 
        # Inizializzo i parametri della visualizzazione.
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.view_init(axIncl[0], axIncl[1])        # rotazione assi
        if not axis:
            ax.axis('off')
        
        # Se curveLevels e' True eseguo visualizzazione della superficie mediante
        # curve di livello.
        if curveLevels:
            # Settaggio dei livelli.
            tol = 6e-3
            (mn,mx) = (min(self.plotHeights),max(self.plotHeights))
            if mn < -tol:
                lvs = int(lvs/2)
                levels = np.block([np.linspace(mn,-tol,lvs),
                                   np.linspace(tol,mx,lvs)])
            else:
                levels = np.linspace(tol,mx,lvs)

            # Preparazione della visualizzazione.
            ax.tricontour(self.plotTriangulation, self.plotHeights, levels=levels, cmap=cmap, linewidths=0.5)

        # Caso di visualizzazione classica.
        else:
            ax.plot_trisurf(self.plotTriangulation, self.plotHeights)
        
        # Costruzione della triangolazione sottostante.
        self.domain.domainPlot(labels=False, legend=False, showExecute=False, color='#000000', colorPS='#000000', linewidths=0.4)
        
        # Salvo l'immagine, se richiesto.
        if save:
                name = input('Nome del file?\n')
                plt.savefig(name, dpi=600)
        plt.show()

        # Mostro l'immagine a schermo.
        plt.show()

    def splineContourPlot(self, cmap='twilight_shifted', save=False, axis=False, showLift=False, lvs=40):
        # Visualizza curve di livello su grafico 2D della superficie B-spline.
        # showLift: se true mostra i triangoli di sollevamento utilizzati nella costruzione della superficie;
        # lvs: numero di livelli da mostrare.
        if not axis:
            plt.axis('off')
        plt.axis('equal')

        # Settaggio dei livelli.
        (mn,mx) = (min(self.plotHeights),max(self.plotHeights))
        tol = 6e-3
        if mn < -tol:
            lvs = int(lvs/2)
            levels = np.block([np.linspace(mn,-tol,lvs),np.linspace(tol,mx,lvs)])
        else:
            levels = np.linspace(tol,mx,lvs)
        
        # Visualizzazione delle curve di livello.
        plt.tricontour(self.plotTriangulation, self.plotHeights, levels=levels, cmap=cmap, linewidths=0.4)
        
        # Visualizzazione del dominio e dei triangoli di sollevamento, se richiesti.
        self.domain.domainPlot(labels=False, legend=False, showExecute=False, color='#000000', colorPS='#000000')
        if showLift:
            self.lift_triangles.liftPlot(onlyLifts=True, legend=False, labels=False, showExecute=False)
        
        # Salvo l'immagine, se richiesto.
        if save:
                name = input('Nome del file?\n')
                plt.savefig(name, dpi=600)
        plt.show()

