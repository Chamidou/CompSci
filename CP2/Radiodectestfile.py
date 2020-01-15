"""
Test file of Radiodec.py
"""
from Radiodec import Radecay


def main():
    nsi = int(input("Please give the array size "))
    deltat = float(input("Please give the time step value "))
    lamb = float(input("Please give the value of the decay constant "))
    proba = deltat*lamb
    matrixtest = Radecay(nsi, deltat, lamb)
    matrixtest.array(deltat, lamb)


main()
