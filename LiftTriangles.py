import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.optimize import LinearConstraint, minimize
from numpy.linalg import norm, solve
from DomainTriangulation import DomainTriangulation

class LiftTriangles:
    def __init__(self, startTr: Triangulation):
        # Inizializzazione delle triangolazioni sul dominio.
        self.domain = DomainTriangulation(startTr)

        # Creazione dei triangoli di sollevamento.
        self.makeLiftTriangles()

    def makeLiftTriangles(self):
        # Inizializzazione liste dei triangoli di sollevamento e relativi vertici.
        VSollTrList = []
        VSollVxList = []
        TSollTrList = []
        TSollVxList = []
        TSollLabels = []       # lista di pedici dei vertici per i triangoli TSoll.
        ESollVxList = []
        ESollLabels = []

        # TSollFinder e' un dizionario tale da associare ad ogni coppia (n,i), dove n e 'un triangolo della triangolazione di partenza e V_i un suo vertice, la posizione in TSollTr del triangolo di sollevamento corrispondente.
        self.TSollFinder = {}

        # Inizializzazione liste dei coefficienti sigma.
        # TSollSigma[m] = sigma_m^t;
        # BESollSigma[i,j] = sigma_(i,j)^e
        self.TSollSigma = {}
        self.ESollSigma = {}

        # Inizializzazione del dizionario dei coefficienti alpha.
        # alpha[i][j] = alpha_i,(j+1) j=0,1,2;
        # alpha[i][j] = alpha^x_i,(j-2) j=3,4,5;
        # alpha[i][j] = alpha^y_i,(j-5) j=6,7,8.
        self.alpha = {}

        # Triangoli di sollevamento associati ai vertici.
        for i in range(len(self.domain.startPoints)):
            # Ricavo il triangolo per il vertice corrente.
            liftTr = self.findVertexLiftTriangle(i)

            # Aggiungo punti e triangolo in lista.
            codTr = []
            for pt in liftTr:
                VSollVxList.append(pt)
                codTr.append(len(VSollVxList)-1)
            VSollTrList.append(np.array(codTr))

        # Creo la triangolazione contenente i triangoli di sollevamento associati ai vertici.
        self.VSollVx = np.array(VSollVxList)
        self.VSollTr = np.array(VSollTrList)
        self.VSoll = Triangulation(self.VSollVx[:,0], self.VSollVx[:,1], self.VSollTr)

        # Triangoli di sollevamento associati alla triangolazione.
        for n in range(len(self.domain.startTr.triangles)):
            for i in self.domain.startTr.triangles[n]:
                codTr = []

                # Aggiungo il punto S_{i,n}^t.
                TSollVxList.append(2/3*self.domain.startPoints[i]+1/3*self.domain.incenters[n])
                TSollLabels.append('$S^t_{'+str(i+1)+','+str(n+1)+'}$')
                codTr.append(len(TSollVxList)-1)

                # Ricerca dei triangoli PS m e m' a vertice V_i.
                mList, jList = self.getPSTriangles(n,i)
                
                # Calcolo i coefficienti sigma^t.
                mx = max(self.domain.xi[n,jList[0]]/self.domain.lamb[jList[0],i] + self.domain.xi[n,jList[1]]/self.domain.lamb[jList[1],i], 1)
                for l in range(2):
                    self.TSollSigma[mList[l]] = self.domain.lamb[jList[l],i]/2 * mx

                # Aggiunta dei punti Q_m^t e Q_m'^t.
                for l in range(2):
                    TSollVxList.append(TSollVxList[len(TSollVxList)-1-l] + 2/3*self.TSollSigma[mList[l]]*(self.domain.startPoints[jList[l]]-self.domain.startPoints[i]))
                    TSollLabels.append('$Q^t_{'+str(mList[l]+1)+'}$')
                    codTr.append(len(TSollVxList)-1)
                
                # Aggiunta del triangolo di sollevamento e aggiornamento di TSollFinder.
                TSollTrList.append(np.array(codTr))
                self.TSollFinder[n,i] = len(TSollTrList)-1

        # Creo la triangolazione contenente i triangoli di sollevamento associati ai vertici.
        self.TSollVx = np.array(TSollVxList)
        self.TSollTr = np.array(TSollTrList)
        self.TSollVxLabels = TSollLabels
        self.TSoll = Triangulation(self.TSollVx[:,0], self.TSollVx[:,1], self.TSollTr)

        # Segmenti di sollevamento associati ai lati di bordo.
        # ESollPos[i,j] memorizza la posizione di S^e_{ij} in ESollVx.
        self.ESollPos = {}
        for edge in self.domain.splitPoints:
            # se il lato edge sta nel bordo si ricavano i punti S^e_ij e Q_{ij}^e 
            if self.domain.splitPoints[edge][1]:
                for i in range(2):
                    # punti S^e_{ij}
                    ESollVxList.append(2/3*self.domain.startPoints[edge[i]] + 1/3*self.domain.splitPoints[edge][0])
                    ESollLabels.append('$S^e_{'+str(edge[i]+1)+','+str(edge[(i+1)%2]+1)+'}$')

                    # aggiorno ESollPos
                    self.ESollPos[edge[i],edge[(i+1)%2]] = len(ESollVxList)-1

                    # coefficienti sigma^e_{ij}
                    self.ESollSigma[edge[i],edge[(i+1)%2]] = self.domain.lamb[edge[(i+1)%2],edge[i]]/2

                    # punti Q^e_{ij}
                    ESollVxList.append(ESollVxList[-1] + 2/3 * self.ESollSigma[edge[i],edge[(i+1)%2]] * (self.domain.startPoints[edge[(i+1)%2]] - self.domain.startPoints[edge[i]]))
                    ESollLabels.append('$Q^e_{'+str(edge[i]+1)+','+str(edge[(i+1)%2]+1)+'}$')

        # Salvo vertici ed etichette dei segmenti di sollevamento al bordo.
        self.ESollVxLabels = ESollLabels
        self.ESollVx = np.array(ESollVxList)

    def findVertexLiftTriangle(self, i: int):
        # Inizializzazione della lista di vincoli LinearConstraint.
        ConstrList = []

        # Vincoli delle coordinate baricentriche.
        A = np.array([[1,1,1,0,0,0,0,0,0],
                      [0,0,0,1,1,1,0,0,0],
                      [0,0,0,0,0,0,1,1,1]])
        lb = ub = np.array([1,0,0])
        ConstrList.append(LinearConstraint(A,lb,ub))

        # Vincolo di non-negativita' degli alfa_i,r.
        A = np.eye(3,9)
        lb = np.zeros(3)
        ub = np.infty * np.ones(3)
        ConstrList.append(LinearConstraint(A,lb,ub))

        # Vincoli associati ai triangoli.
        for tr_index,_,_ in self.domain.TList[i]:
            tr = self.domain.PS_triangles[tr_index]

            # Costruzione della matrice di vincolo LinearConstraint.
            def matr(diffr, diffz):
                return np.block([[np.eye(3), diffr[0]/3 * np.eye(3), diffr[1]/3 * np.eye(3)],
                                 [np.eye(3), diffz[0]/3 * np.eye(3), diffz[1]/3 * np.eye(3)]])
                
            # Creazione del vincolo associato a tr e inserimento in lista.
            lb = np.zeros(6)
            ub = np.infty * np.ones(6)
            Vi = self.domain.PS_points[tr[0]]
            Zn = self.domain.PS_points[tr[1]]
            Rij = self.domain.PS_points[tr[2]]
            ConstrList.append(LinearConstraint(matr(Rij-Vi,Zn-Vi),lb,ub))

        # vettore iniziale
        a0 = np.array([1/3,1/3,1/3,1,-1,0,2,-1,1])

        # funzione obiettivo
        def f(a):
            return -a[3]*a[7]+a[6]*a[4]
        
        # Esecuzione del problema di ottimizzazione per la ricerca del triangolo di area minima.
        alpha = minimize(f,a0,constraints=ConstrList).x

        # Memorizzazione dei coefficienti alpha.
        self.alpha[i] = alpha
        
        # Ricavo le coordinate del triangolo di sollevamento risolvendo opportuno sistema lineare.
        M = np.block([[alpha[0:3]],[alpha[3:6]],[alpha[6:9]]])
        B = np.block([[self.domain.startPoints[i], 1],[1,0,0],[0,1,0]])
        Q = solve(M,B)

        # Restituisco le coordinate Qi,1 Qi,2 Qi,3 ottenute.
        return (Q[0,0:2], Q[1,0:2], Q[2,0:2])

    def getPSTriangles(self, n: int, i: int):
        # Fornisce (mList, jList): mList e' la lista delle posizioni
        # in PS_Triangles dei triangoli di PS aventi Vi come vertice, jList e' la lista 
        # dei Vj adiacenti.
        mList = []
        jList = []
        for el in self.domain.TList[i]:
                if norm(el[2]-self.domain.incenters[n]) == 0:
                    mList.append(el[0])
                    jList.append(el[1])
        return (mList, jList)

    def liftPlot(self, showTitle=False, labels=True, legend=True, showExecute=True, fontsize=7, \
         save=False, onlyLifts=False):
        # Visualizza i triangoli di sollevamento sul dominio.
        # showTitle: se True mostra il titolo;
        # labels: se True mostra nomi dei punti;
        # legend: se True mostra la legenda;
        # fontsize: la dimensione del font per i punti;
        # showExecute: se True esegue plt.show() dentro il metodo;
        # fontsize: la dimensione del font per i punti;
        # save: se True salva l'immagine;
        #onlyLifts: se True mostra soltanto le strutture di sollevamento.
        if showTitle:
            plt.title('Triangoli di sollevamento per la costruzione di una base $\mathbb{S}_3^1(\Delta_{PS})$.')
        
        # Se richiesta, mostro triangolazione sottostante.
        if not onlyLifts:
            self.domain.domainPlot(legend=False, showExecute=False, labels=False, color='#c6c5c6', \
            colorPS='#c6c5c6', showDP=True)

        # aggiustamento degli assi
        plt.axis('equal')

        # triangoli di sollevamento relativi ai vertici
        plt.triplot(self.VSoll, color='#fd3612', linewidth=0.9, label='triangoli $\mathcal{Q}_i^v$')

        # triangoli di sollevamento relativi ai traingoli PS
        plt.triplot(self.TSoll, color='#0e03fb', linewidth=0.9, label='triangoli $\mathcal{Q}_{i,n}^t$')

        # segmenti di sollevamento ai lati di bordo
        for i in range(int(len(self.ESollVx)/2)):
            plt.plot(self.ESollVx[2*i:2*(i+1),0], self.ESollVx[2*i:2*(i+1),1], \
            label='segmenti $\mathcal{Q}_{i,j}^e$', color='#fac900', linewidth=0.9)

        # Etichettatura dei punti, se richiesta.
        if labels:
            incr= 0.01
            # vertici
            for v in range(len(self.domain.startPoints)):
                x = self.domain.startPoints[v,0]
                y = self.domain.startPoints[v,1]
                plt.text(x+incr,y+incr,'$V_{'+str(v+1)+'}$',fontsize=fontsize)

            # punti associati ai vertici
            for i in range(len(self.VSollTr)):
                for j in range(len(self.VSollTr[i])):
                    x = self.VSollVx[self.VSollTr[i,j]][0]
                    y = self.VSollVx[self.VSollTr[i,j]][1]

                    # non visualizzo etichetta quando ho coincidenza con un vertice
                    if norm(self.VSollVx[self.VSollTr[i,j]]-self.domain.startPoints[i]) > 1e-2:
                        plt.text(x+incr,y+incr,'$Q^v_{'+str(i+1)+','+str(j+1)+'}$',fontsize=fontsize)

            # punti associati ai triangoli PS
            for i in range(len(self.TSollVxLabels)):
                x = self.TSollVx[i][0]
                y = self.TSollVx[i][1]
                plt.text(x+incr,y+incr,self.TSollVxLabels[i],fontsize=fontsize)

            # segmenti al bordo
            for i in range(len(self.ESollVxLabels)):
                x = self.ESollVx[i][0]
                y = self.ESollVx[i][1]
                plt.text(x-4*incr,y-4*incr,self.ESollVxLabels[i],fontsize=fontsize)

        # Esecuzioni di legenda, plt.show() e salvataggio, se richieste.
        if showExecute:
            if legend:
                (h,l) = plt.gca().get_legend_handles_labels()
                plt.legend(handles=[h[5],h[7],h[9]],labels=[l[5],l[7],l[9]], fontsize='x-small')
            if save:
                name = input('Nome del file?\n')
                plt.savefig(name, dpi=600)
            plt.show()
