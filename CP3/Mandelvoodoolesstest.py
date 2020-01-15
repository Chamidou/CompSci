"""
Test file for Mandelbrot
"""
from Mandelvoodooless import Mandelbrot


def main():
    lim = int(input("Enter the maximum amount of iterations "))
    creamin = -2.025        #Minimum limit for the real part of c
    creamax = 0.6       #Maximum limit for the real part of c
    cimamin = -1.125        #Minium limit for the imaginary part of c
    cimamax = 1.125         #Maximum limit for the imaginary part of c
    mandeltest = Mandelbrot(lim, creamin, creamax, cimamin, cimamax)
    mandeltest.mandeloop()

main()
