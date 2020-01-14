"""
Test file for the project
"""

from Classfile import Solarsys

def main():
    itera = 40000
    deltat = 5000
    sim = Solarsys(itera, deltat)
    sim.readfile(), sim.append(), sim.sumacceleration(), sim.kinetic(), sim.orbitalperiod()
    sim.closeapproachmars()
    sim.closeapproachjupiter()
    sim.display()
main()
