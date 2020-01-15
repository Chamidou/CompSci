'''
Douglas Coenen
Computer Simulation: Checkpoint 2: Radioactive decay
'''
import matplotlib.pyplot as plt
import random
import numpy as np
import math

class Radecay(object):
    def __init__(self, nsi, deltat, lamb):	 #Creates the constructor
         self.length = nsi
         self.atarray = np.zeros((nsi, nsi))        #Creates the matrix full of zeros
         self.proba = deltat * lamb
    def array(self, deltat, lamb):

        count = 0	#Set up a counter for all the loops made
        decayed = 0     #Set up a count for each time an atom decays
        halfdec = (self.length**2)/2        #Sets the amount of atom to stop the loop at half the total
        while decayed < halfdec:
            count += 1      #Every loop, the count holder gains 1
            for i in range(0, self.length):     #Goes through the first dimension of the matrixtest
                for j in range(0, self.length):     #Goes through the lists inside the motherlist
                    if self.atarray[i, j] == 0:     #If the atom is undecayed apply next
                        if random.random() <= self.proba:	#If the generated number is inferior to the calculated proba
                            self.atarray[i, j] = 1	#change the state of the atom into decayed
                            decayed += 1        #add on to the decayed atoms holder
        simhalflife = count * deltat	   #Calculates the simulated halflife
        totalatom = self.length**2
        undec = totalatom - decayed
        theohalflife = (math.log(2))/lamb	#Calculate the theoretical half-life with the given constants
        print("Simulated halflife : " + str(simhalflife) + " min")
        print("Theoretical halflife : " + str(theohalflife) + " min")
        print("Total number of atoms " + str(self.length**2))
        print("Number of undecayed atoms " + str(undec))
        plt.matshow(self.atarray, cmap="hot")       #matshow was used the array representation
        plt.xlabel("Black = Undecayed, White = Decayed ")
        plt.title("Visual represention of decayed and undecayed atoms")
        plt.show()
