import numpy as np
import matplotlib.pyplot as plt 
import random as rnd
from matplotlib.tri import Triangulation
from numpy.linalg import norm

class DomainTriangulation:
    def __init__(self, startTr: Triangulation):
        self.startTr = startTr

        # Creazione della triangolazione raffinata.
        self.PowellSabin()

    def PowellSabin(self):
        # Esegue raffinamento di Powell-Sabin.
        # lamb, mu, xi sono dizionari che memorizzano le coordinate baricentriche necessarie alla costruzione della Powell-Sabin B-spline.
        # R_ij = lamb[i,j] * V_i + lamb[j,i] * V_j;
        # R_ij = mu[n,n'] * Z_n + mu[n',n] * Z_n';
        # Z_n = xi[n,i] * V_i + xi[n,j] * V_j + xi[n,k] * V_k.
        self.lamb = {}
        self.mu = {}
        self.xi = {}
        
        # incenters e splitPoints memorizzano rispettivamente gli incentri e punti di suddivisione sui lati.
        # incenters[n] = Z_n;
        # splitPoints[i,j] = (R_ij, b) con i < j e b==True <=> V_iV_j e' un lato di bordo;
        self.incenters = {}
        self.splitPoints = {}

        # TList[i] e' una lista di elementi della forma (tpos, j, in): tpos e' la posizione in PS_triangles di un triangolo di Powell-Sabin avente V_i come vertice, j e' l'indice del corrispettivo R_ij e in l'incentro Z_n.
        self.TList = {}

        # Salvo i punti della triangolazione iniziale in un solo array.
        self.startPoints = np.array(list(zip(self.startTr.x,self.startTr.y)))

        # startEges e' un dizionario che memorizza per ogni triangolo la lista di lati (i,j) con i<j tali che il lato V_iV_j stia nel triangolo n.
        self.startEdges = {}

        for n in range(len(self.startTr.triangles)):
            # Salvo le coordinate baricentriche dell'incentro qualora non abbia ancora trattato il presente triangolo.
            self.saveIncenter(n)

            # Salvo le coordinate baricentriche dei punti di suddivisione.
            self.saveSplitPoints(n)

        # Inizializzo lista di triangoli e punti del raffinamento PS.
        PS_triangles = []
        trRow = 0
        PS_points = []
        ptRow = 0

        # Colloco i nuovi triangoli: il dizionario pointsPos memorizza le posizioni-riga (ptRow) in PS_points dei punti gia' collocati.
        pointsPos = {}
        for n in range(len(self.startTr.triangles)):
            # collocamento dell'incentro
            PS_points.append(self.incenters[n])
            pointsPos['in' + str(n)] = ptRow
            ptRow += 1

            for edge in self.startEdges[n]:
                # collocamento del punto di suddivisione
                if not edge in pointsPos:
                    PS_points.append(self.splitPoints[edge][0])
                    pointsPos['ed' + str(edge)] = ptRow
                    ptRow += 1

                for v in edge:
                    # collocamento dei vecchi vertici
                    if not v in pointsPos:
                        PS_points.append(self.startPoints[v])
                        pointsPos['pt' + str(v)] = ptRow
                        ptRow += 1 

                    # collocamento dei nuovi triangoli e aggiornamento TList.
                    PS_triangles.append(np.array([pointsPos['pt'+str(v)],pointsPos['in'+str(n)],pointsPos['ed'+str(edge)]]))
                    if not v in self.TList:
                        self.TList[v] = []
                    self.TList[v].append((trRow, edge[(edge.index(v)+1)%2], self.incenters[n]))
                    trRow += 1

        # Creo l'oggetto triangolazione PS-raffinata.
        self.PS_points = np.array(PS_points)
        self.PS_triangles = np.array(PS_triangles)
        self.PowellSabinTr = Triangulation(np.array(PS_points)[:,0],np.array(PS_points)[:,1],np.array(PS_triangles))

    def saveIncenter(self, n: int):
        # Esegue l'inserimento dell'incentro per il triangolo n-esimo solo
        # se non e' gia' stato effettuato.
        if not (n,self.startTr.triangles[n,0]) in self.xi:
            # Vertici del triangolo corrente.
            v0 = self.startPoints[self.startTr.triangles[n,0]]
            v1 = self.startPoints[self.startTr.triangles[n,1]]
            v2 = self.startPoints[self.startTr.triangles[n,2]]

            # Coordinate baricentriche dell'incentro.
            perim = norm(v1-v0) + norm(v2-v1) + norm(v0-v2)
            self.xi[n,self.startTr.triangles[n,0]] = norm(v2-v1)/perim
            self.xi[n,self.startTr.triangles[n,1]] = norm(v0-v2)/perim
            self.xi[n,self.startTr.triangles[n,2]] = norm(v1-v0)/perim

            # Coordinate cartesiane dell'incentro.
            self.incenters[n] = self.xi[n,self.startTr.triangles[n,0]]*v0 + self.xi[n,self.startTr.triangles[n,1]]*v1 + self.xi[n,self.startTr.triangles[n,2]]*v2

    def saveSplitPoints(self, n: int):
        # Esegue l'inserimento dei punti di suddivisione per il triangolo n-esimo.
        for tr_index in range(3):
            # Vertici del lato corrente.
            i = min(self.startTr.triangles[n,tr_index],self.startTr.triangles[n,(tr_index+1)%3])
            j = max(self.startTr.triangles[n,tr_index],self.startTr.triangles[n,(tr_index+1)%3])
            vi = self.startPoints[i]
            vj = self.startPoints[j]

            # Aggiorno startEdges.
            if not n in self.startEdges:
                self.startEdges[n] = []
            self.startEdges[n].append((i,j))

            # Lato al bordo.
            if self.startTr.neighbors[n,tr_index] == -1:
                self.lamb[i,j] = 1/2
                self.lamb[j,i] = 1/2
                self.splitPoints[i,j] = (1/2 * vi + 1/2 * vj, True)
            
            # Lato interno.
            else:
                # Ricerca del triangolo adiacente ed inserimento del relativo incentro, se necessario.
                n_prime = self.startTr.neighbors[n,tr_index]
                self.saveIncenter(n_prime)

                # Selezione degli incentri coinvolti nella suddivisione.
                zn = self.incenters[n]
                zn_prime = self.incenters[n_prime]

                # Calcolo delle coordinate baricentriche del punto di suddivisione.
                diffz = zn_prime - zn
                diffv = vi - vj

                if diffz[0] != 0:
                    self.lamb[j,i] = (vi[1]*diffz[0]-vi[0]*diffz[1]+zn[0]*zn_prime[1]- zn_prime[0]*zn[1])/(diffz[0]*diffv[1]-diffz[1]*diffv[0])
                    self.lamb[i,j] = 1 - self.lamb[j,i]
                    self.mu[n_prime,n] = (vi[0]-zn[0]-self.lamb[j,i]*diffv[0])/diffz[0]
                    self.mu[n,n_prime] = 1 - self.mu[n_prime,n]
                else:
                    self.lamb[j,i] = (vi[0]*diffz[1]-vi[1]*diffz[0]+zn[1]*zn_prime[0]- zn_prime[1]*zn[0])/(diffz[1]*diffv[0]-diffz[0]*diffv[1])
                    self.lamb[i,j] = 1 - self.lamb[j,i]
                    self.mu[n_prime,n] = (vi[1]-zn[1]-self.lamb[j,i]*diffv[1])/diffz[1]
                    self.mu[n,n_prime] = 1 - self.mu[n_prime,n]

                # Calcolo del punto di suddivisione.
                self.splitPoints[i,j] = (self.lamb[i,j] * vi + self.lamb[j,i] * vj, False)

    def getMultiIndices(self):
        # Metodo che genera tutti i multi-indici per spline cubiche.
        return [(i,j,k) for i in range(4) for j in range(4) for k in range(4) if i+j+k==3]

    def getZnOpposites(self, n:int):
        # Metodo che restiuisce una lista di tuple (m,m'), dove m e m' sono posizioni in PS_triangles di triangoli che condividono un lato ZnRij, e una lista di tuple (i,j), ciascuna corrispondente all'Rij condiviso.
        res1 = []
        res2 = []
        for i in self.startTr.triangles[n]:
            for el in self.TList[i]:
                if norm(el[2]-self.incenters[n]) == 0:
                    res1.append((el[0],self.PowellSabinTr.neighbors[el[0],1]))
                    res2.append((i,el[1]))
        return (res1, res2)

    def getRijOpposite(self, n: int, edge: tuple):
        # Metodo che restituisce la posizione m' del triangolo PS giacente sul lato edge=(i,j), avente in comume il vertice Vi con il triangolo PS di partenza.
        n_prime = self.getCommonTriangle(n, edge)
        (i,j) = edge
        for el in self.TList[i]:
            if el[1] == j and norm(el[2]-self.incenters[n_prime]) == 0:
                return el[0]
        raise ValueError('Edge lato di bordo.')
        
    def getCommonTriangle(self, n: int, edge: tuple):
        # Metodo che restituisce la posizione in startTr.triangles del triangolo avente in comune il lato edge=(i,j) con il triangolo in posizione n dato.
        (i,j) = edge
        if not i in self.startTr.triangles[n] or not j in self.startTr.triangles[n]:
            raise ValueError('Inserted edge doesn not belong to triangle.')

        tri = self.startTr.triangles[n]
        (ind1,ind2) = (min(list(tri).index(i),list(tri).index(j)), max(list(tri).index(i),list(tri).index(j)))
        if (ind1,ind2) == (0,2):
            return self.startTr.neighbors[n,2]
        else:
            return self.startTr.neighbors[n,ind1]

    def domainPlot(self, showTitle=False, showPS=True, showDP=False, labels=True, legend=True, showExecute=True, linewidths=0.6, fontsize=8, color='#5a5e99', colorPS='#969dff', save=False):
        # Visualizza il dominio della spline con relative triangolazioni.
        # showTitle: se True mostra il titolo;
        # showPS: se True mostra raffinamento di Powell-Sabin;
        # showDP: se True mostra proiezioni sul dominio dei punti di controllo;
        # labels: se True mostra nomi dei punti;
        # legend: se True mostra la legenda;
        # showExecute: se True esegue plt.show() dentro il metodo;
        # linewidths: spessore delle linee;
        # fontsize: la dimensione del font per i punti;
        # color: il colore per la triangolazione originale;
        # colorPS: il colore per la triangolazione Powell-Sabin;
        # save: se True, salva l'immagine.
        if showTitle:
            plt.title('Triangolazione iniziale e suddivisione Powell-Sabin.')
        
        if showExecute:
            plt.axis('equal')
        plt.axis('off')

        # Visualizzazione della triangolazione originale e, se richiesto,
        # del relativo raffinamento.
        plt.triplot(self.startTr, color=color, linewidth=linewidths, label='triangolazione $\Delta$')
        if showPS:
            plt.triplot(self.PowellSabinTr, linestyle='dotted', linewidth=linewidths, color=colorPS, label='triangolazione $\Delta_{PS}$')

        # Etichettatura dei punti, se richiesta.
        if labels:
            incr= 0.01
            # incentri
            for n in self.incenters:
                x = self.incenters[n][0]
                y = self.incenters[n][1]
                plt.text(x+incr,y+incr,'$Z_{'+str(n+1)+'}$',fontsize=fontsize)

            # vertici
            for v in range(len(self.startPoints)):
                x = self.startPoints[v,0]
                y = self.startPoints[v,1]
                plt.text(x+incr,y+incr,'$V_{'+str(v+1)+'}$',fontsize=fontsize)

            # punti suddivisione
            for r in self.splitPoints:
                x = self.splitPoints[r][0][0]
                y = self.splitPoints[r][0][1]
                plt.text(x+incr,y+incr,'$R_{'+str(r[0]+1)+','+str(r[1]+1)+'}$',fontsize=fontsize)

        # Visualizzazione di punti di controllo (con IDM=insieme di determinazione minimo),
        # se richiesta.
        if showDP:
            multiInd = self.getMultiIndices()
            
            # Inizializzazione liste di punti di controllo e insieme di determinazione minimo.
            DP = []
            MDS = []

            for v in range(len(self.startPoints)):
                # Per ogni vertice un triangolo adiacente e' scelto a caso.
                trChosen = self.TList[v][rnd.randint(0,len(self.TList[v])-1)][0]

                for tr_index, w, _ in self.TList[v]:
                    triangle = self.PS_triangles[tr_index]
                    # Veritci del triangolo corrente.
                    p1 = self.PS_points[triangle[0]]
                    p2 = self.PS_points[triangle[1]]
                    p3 = self.PS_points[triangle[2]]
                    
                    for (i,j,k) in multiInd:
                        # i vertici della triangolazione originale sono nell'IDM
                        if (i,j,k) == (3,0,0):
                            MDS.append(p1)

                        # i baricentri dei triangoli PS sono nell'IDM
                        elif (i,j,k) == (1,1,1):
                            MDS.append(1/3*p1 + 1/3*p2 + 1/3*p3)

                        # i punti a indici (2,1,0) e (2,0,1) per un triangolo scelto a caso sono nell'IDM
                        elif tr_index == trChosen and ((i,j,k) == (2,1,0) or (i,j,k) == (2,0,1)):
                            MDS.append(i/3*p1 + j/3*p2 + k/3*p3)

                        # i punti a indici (1,0,2) se appartenenti ad un lato di bordo sono nell'IDM
                        elif self.splitPoints[min(v,w),max(v,w)][1] and (i,j,k) == (1,0,2):
                            MDS.append(i/3*p1 + j/3*p2 + k/3*p3)

                        # tutti gli altri sono normali punti di controllo
                        else:
                            DP.append(i/3*p1 + j/3*p2 + k/3*p3)

            # Salvataggio dei punti di controllo e dell'IDM.
            self.DP = np.array(DP)
            self.MDS = np.array(MDS)

            # Visualizzazione dei punti di controllo e IDM.
            plt.plot(self.DP[:,0],self.DP[:,1], color=colorPS, linestyle='', marker='o', markersize=2)
            plt.plot(self.MDS[:,0],self.MDS[:,1], color=color, linestyle='', marker='o', markersize=2, label='ins. determ. minimo')
                        
        # Esecuzioni di legenda, plt.show() e salvataggio, se richieste.
        if showExecute:
            if legend:
                (h,l) = plt.gca().get_legend_handles_labels()
                if showDP:
                    plt.legend(handles=[h[0],h[2],h[4]],labels=[l[0],l[2],l[4]], fontsize='x-small')
                else:
                    plt.legend(handles=[h[0],h[2]],labels=[l[0],l[2]], fontsize='x-small')
            if save:
                name = input('Nome del file?\n')
                plt.savefig(name, dpi=600)
            plt.show()

                
            

