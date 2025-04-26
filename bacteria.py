from fastaReader import fastaReader
import random
import numpy 
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria():
    
    
    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction =0
        self.NFE = 0
        
    def showGenome(self):
     for seq in self.matrix.seqs:
        print(seq)

    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = numpy.array(copy.deepcopy(self.matrix.seqs))
        return newBacteria

    def tumboNado(self, numGaps, tablaFitness=None, mejorHistorico=None): # Actualizado
        # Se agrego influencia historica y la mutación adaptativa
        self.cuadra()
        matrixCopy = copy.deepcopy(self.matrix.seqs)
        matrixCopy = matrixCopy.tolist()
        num_secuencias = len(matrixCopy)
        longitud_alineamiento = len(matrixCopy[0]) if matrixCopy else 0
        prob_influencia_historica = 0.1

        if tablaFitness is None:
            gapRandomNumber = random.randint(0, numGaps)
            for _ in range(gapRandomNumber):
                seqnum = random.randint(0, num_secuencias - 1)
                pos = random.randint(0, longitud_alineamiento)
                matrixCopy[seqnum] = matrixCopy[seqnum][:pos] + "-" + matrixCopy[seqnum][pos:]
        else:
            max_fitness = max(tablaFitness) if tablaFitness else 1
            min_fitness = min(tablaFitness) if tablaFitness else 0
            rango_fitness = max_fitness - min_fitness if max_fitness > min_fitness else 1
            fitness_relativo = (self.fitness - min_fitness) / rango_fitness if rango_fitness > 0 else 0.5
            prob_mutacion = 1 - fitness_relativo

            num_gaps_adaptativo = max(1, int(numGaps * (1 - fitness_relativo)))
            for _ in range(num_gaps_adaptativo):
                if mejorHistorico is not None and random.random() < prob_influencia_historica and longitud_alineamiento > 0:
                    columna_aleatoria = random.randint(0, longitud_alineamiento - 1)
                    if all(seq[columna_aleatoria] == '-' for seq in mejorHistorico):
                        secuencia_aleatoria = random.randint(0, num_secuencias - 1)
                        posicion_insercion = columna_aleatoria
                        matrixCopy[secuencia_aleatoria] = matrixCopy[secuencia_aleatoria][:posicion_insercion] + "-" + matrixCopy[secuencia_aleatoria][posicion_insercion:]
                    else:
                        seqnum = random.randint(0, num_secuencias - 1)
                        pos = random.randint(0, longitud_alineamiento)
                        matrixCopy[seqnum] = matrixCopy[seqnum][:pos] + "-" + matrixCopy[seqnum][pos:]
                else:
                    seqnum = random.randint(0, num_secuencias - 1)
                    pos = random.randint(0, longitud_alineamiento)
                    matrixCopy[seqnum] = matrixCopy[seqnum][:pos] + "-" + matrixCopy[seqnum][pos:]

        self.matrix.seqs = numpy.array(matrixCopy)
        self.cuadra()
        self.limpiaColumnas()        


    def cuadra(self):
        # Se rellena con gaps las secuencias mas cortas.
        import numpy
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-"*(maxLen-len(seq[i]))
        self.matrix.seqs = numpy.array(seq)
        

    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True
    


    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1
        
            
    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]





    def getColumn(self, col):
        column = []
        for i in range(len(self.matrix.seqs)):
            column.append(self.matrix.seqs[i][col])
        return column
    


    def autoEvalua(self):   
        evaluador = evaluadorBlosum()
        score = 0
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            gapCount = column.count("-")
            column = [x for x in column if x != "-"]
            pares = self.obtener_pares_unicos(column)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])
            score -= gapCount*2
        self.blosumScore = score
        self.NFE += 1
        

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)   