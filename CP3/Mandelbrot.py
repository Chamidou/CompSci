
"""
Douglas Coenen
Mandelbrot set
"""

import matplotlib.pyplot as plt
import numpy as np


class Mandelbrot(object):

    def __init__(self, lim, creamin, creamax, cimamin, cimamax):
        self.lim = lim
        self.c_reamin = creamin
        self.c_reamax = creamax
        self.c_imamin = cimamin
        self.c_imamax = cimamax
        v = self.lim * 1j       #Multiply the limit by 1j to create a complex number
        x,y=np.ogrid[creamin:creamax:v , cimamin:cimamax:v]     #creates meshgrid of set limits and the number of steps is the input limit
        self.c = x + 1j*y       #Puts "c" as an array of x as the real part and y the imaginary part

    def mandeloop(self):

        z = 0       #Define the starting point of z

        for i in range(255):
            z = z**2 + self.c       #Apply the function Z(n+1) = Z**2 + c

        stop = z**2 < 4        #Vectorizes the array to iterate each elements one by one
        plt.imshow(stop.transpose(), extent=[self.c_reamin,self.c_reamax, self.c_imamin,self.c_imamax], cmap="hot")
        plt.show()
