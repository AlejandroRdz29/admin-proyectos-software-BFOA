from bacteria import bacteria
from quimiotaxis import quimiotaxis

import numpy


def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction
    
def validaSecuencias(path, veryBest):
    #clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    #descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-","")
    #tempBacteria.tumboNado(1)    

    #valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

if __name__ == "__main__":
    poblacion = []
    path = "multifasta.fasta"
    numeroDeBacterias = 6
    numRandomBacteria = 2
    iteraciones = 5
    tumbo = 3                                            
    nado = 4
    chemio = quimiotaxis()
    veryBest = bacteria(path)  
    veryBest_historico = None # Nuevo parametro                      
    tempBacteria = bacteria(path)                      
    original = bacteria(path)                           
    globalNFE = 0                                      

    dAttr= 0.15     #0.1 
    wAttr= 0.003    #0.2 
    hRep=dAttr
    wRep= 0.002    #10

    for i in range(numeroDeBacterias):                                       
        poblacion.append(bacteria(path))


    for iteracion in range(iteraciones):
        factor = iteracion / iteraciones 
        
        # Adaptación de parámetros dinamicamente (Nuevo)
        dAttr_iter = dAttr * (1 - factor)  
        wAttr_iter = wAttr * (1 - factor)
        hRep_iter = hRep * (1 + factor)    
        wRep_iter = wRep * (1 + factor)

        tabla_fitness = [bact.fitness for bact in poblacion]

        # Tumbo-Nado por cada bacteria con mutación adaptativa e influencia histórica (Nuevo)
        for bact in poblacion:
            if veryBest_historico is not None:
                bact.tumboNado(tumbo, tabla_fitness, veryBest_historico.matrix.seqs) # Pasa el mejor histórico
            else:
                bact.tumboNado(tumbo, tabla_fitness)
            bact.autoEvalua()

        # Aplicar quimiotaxis con los parámetros ajustados
        chemio.doChemioTaxis(poblacion, dAttr_iter, wAttr_iter, hRep_iter, wRep_iter)

        # Actualizar el conteo global de evaluaciones de la función objetivo
        globalNFE += chemio.parcialNFE
        best = max(poblacion, key=lambda x: x.fitness)
        
        # Si se encuentra una mejor bacteria, se clona
        if (veryBest == None) or (best.fitness > veryBest.fitness):
            clonaBest(veryBest, best)
            veryBest_historico = veryBest.clonar(path) # Actualiza el mejor histórico (Nuevo)


        print("Iteración: ", iteracion, "Interacción: ", veryBest.interaction, "Fitness: ", veryBest.fitness, " NFE:", globalNFE)
        
        # Ajustar el número de bacterias eliminadas dinámicamente
        eliminar_porcentaje = 0.5 - (0.3 * factor)
        chemio.eliminarClonar(path, poblacion)
        
        # Insertar bacterias aleatorias para mejorar la diversidad
        chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)
        
        # Se imprime la población
        # print("Población: ", len(poblacion))

    # Se imprime la mejor bacteria
    veryBest.showGenome()
    validaSecuencias(path, veryBest)