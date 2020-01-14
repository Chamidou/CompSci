"""
Class file for CompSi project.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm


class Solarsys(object):
    def __init__(self,itera,deltat):    #Sets up the constructor by importing the number of iterations and the number of timesteps
        self.itera = itera      #Make the constants accessible throughout the program by adding slef.
        self.deltat = deltat
        self.G = 6.67408e-11
        self.numlines = 0       #Setting the number of lines counted in the text file to zero, will be used in the next method
        self.bodyname = []      #List which holds the names of the planets

    def readfile(self):
        file_in = open("infofile.txt", "r")     #Import the text file and reads
        data_input = file_in.readlines()        #Read each lines
        data = []           #Creates an empty list for storing the data

        for line in data_input:     #For loop which goes through each line of the text file
            if not line.startswith("#"):       #Adds a recognition of the "#" sign as a comment in the text file and thus ignore the line which starts with it
                tokens = line.split(", ")      #Split the lines into "tokens", so everything between comas will be treated individually
                self.numlines += 1      #Count the lines and in the case, each line is a different planet/satellite. Indicates the total number of bodies
                self.bodyname.append(tokens[0])     #Append the first token of the line to the bodyname list as the first token is the name
                for i in range(1, len(tokens)):     #Creates a for loop which goes through the line and append each token to the list.
                    data.append(float(tokens[i]))
        self.data_array = np.array(data)        #Transform the list into a numpy array
        self.data_array = self.data_array.reshape(self.numlines, 5)     #Reshape the data array so that one line represents one planet with
                                                                        #the 5 parametre (Mass, Position(x,y), Velocity(x,y)

    def append(self):
        self.masslist = np.zeros((self.numlines))        #Create a numpy array which will hold the masses
        self.wholedata = np.zeros((self.itera, self.numlines, 3, 2))    #Create the 4D np array dataset
        for planet in range(self.numlines):
            self.masslist[planet] = self.data_array[planet, 0]      #Append the mass of each body to a separate list
            for j in range(5):
                if j == 1:      #Because the list is composed of x and y coordinates, x component of Position is at index 1 and the y at 2
                    self.wholedata[0, planet, 0] = np.array([self.data_array[planet, j], self.data_array[planet, j + 1]])
                if j == 3:      #So here it starts at 3 since this is where there is the velocity coordinates
                    self.wholedata[0, planet, 1] = np.array([self.data_array[planet, j], self.data_array[planet, j + 1]])
#Here, the planet index go from 0 to 7 as this, Sun, Mercury, Venus, Earth, Mars, Jupiter, SatelliteM, SatelliteJ
#Each have a sub array with three components: Position, Velocity, Acceleration. Index 0 is the position, index 1 is velocity
#and index 2 is acceleration. Each of them have a subarray size 2 which holds their value in x and y
    def sumacceleration(self):      #Main method
        data = self.wholedata       #Refer to self.wholedata as just data for convenience, as it is way shorter to write and will be used extensively
        n = 0       #Set a counter to 0 for the satellite velocity measurement at the end of the method.
        for itera in tqdm(range(self.itera-1)):       #For each iterations: ... ... until the total -1 since Beeman appends at i+1

            if itera == 0:      #For the first iteration:
                for planet in range(self.numlines):     #for each planet:
                    sumaccel = np.array([0.0, 0.0])     #Resets the sum of acceleration to 0
                    for k in range(self.numlines):      #This for loop is designed to use k to go through all the Other planets,
                                                        #so if planet = 1 (Mercury), k = 0,2,3,4,5,6,7
                        if k != planet:     #Refers to the accelcalc method to calculate acceleration the position of two planets and the
                            accele = self.accelcalc(data[itera, planet, 0], data[itera, k, 0], self.masslist[k])    #mass of the other one
                            sumaccel += accele      #Add each acceleration together to get the total

                    sumaccel *= -self.G     #Multiply by -G as seen in the equations in the report
                    data[itera, planet, 2] = sumaccel       #Append this new value to the planet it belongs to

                    #Uses the Euler-Cromer method to calculate the Position and Velocity of the next iteration (itera + 1)
                    #using the present Position, Velocity and acceleration.
                    data[itera + 1, planet, 0], data[itera + 1, planet, 1] = self.euler(planet, data[itera, planet, 0],
                                                                                                    data[itera, planet, 1])

                for planet in range(self.numlines):     #This second loop is for calculating the acceleration of the next step
                    sumaccel = np.array([0.0, 0.0])     #where itera = 1.
                    for k in range(self.numlines):

                        if k != planet:     #This uses the same method than the previous loop, appart from using the value from the next step (i+1)
                            accele = self.accelcalc(data[itera + 1, planet, 0], data[itera + 1, k, 0], self.masslist[k])
                            sumaccel += accele
                    sumaccel *= -self.G
                    data[itera + 1, planet, 2] = sumaccel   #And appends it to i+1 slot for acceleration
                    #So now the initial conditions and the first iteration has all Positions, Velocities and Accelerations calculated.

            if itera > 0:       #For the rest of the simulation, it will use the Beeman scheme.
                for planet in range(self.numlines):     #First loop is designed to go through all planet to calculate their new position
                                                        #at i+1 using velocity and the acceleration from the current and previous timestep
                    data[itera + 1, planet, 0] = data[itera, planet, 0] + data[itera, planet, 1] * self.deltat + ((
                            4 * data[itera, planet, 2] - data[itera - 1, planet, 2]) / 6.0) * (self.deltat ** 2)
                    #Equation used is present in the index of report for a well written version.

                for planet in range(self.numlines):     #This loop goes through all planets to calculate their acceleration
                    sumaccel = np.array([0.0, 0.0])     #from the freshly calculated position verctors
                    for k in range(self.numlines):      #using the same method than the previous step
                        if planet != k:

                            accele = self.accelcalc(data[itera + 1, planet, 0], data[itera + 1, k, 0], self.masslist[k])
                            sumaccel += accele

                    sumaccel *= -self.G
                    data[itera + 1, planet, 2] = sumaccel       #Appends it at the next iteration (i + 1)

                for planet in range(self.numlines):     #This last for loop if for calculating the new velocity. Using the current
                                                        #position and velocity, as well as the acceleration from the next, current and previous iteration
                    data[itera + 1, planet, 1] = data[itera, planet, 1] + (
                            2 * data[itera + 1, planet, 2] + 5 * data[itera, planet, 2] - data[itera - 1, planet, 2]) * self.deltat / 6.0
                #The order of Calculation (Position, Acceleration and Velocity is intended as Velocity + 1 needs Acceleration + 1
                #which in turn needs the Position + 1. But Position + 1 does not need Acceleration + 1 nor Velocity + 1, so it comes first.

                if itera*self.deltat >= 6.3e7 and itera*self.deltat <= 6.35e7: #If the elapsed time is between 2 years and 2 years and 6 days
                #to simulate a constant burn and it is easier to deal with a few extra iterations with low speed than one single with high speed
                   data[itera + 1, 6, 1] = np.array((0.0, 33.111e3))    #Replaces the Velocity vector by a set value for a few iterations

                if itera*self.deltat >= 9.45e7 and itera*self.deltat <= 9.5e7:#Does the same as the above but start 3 years from the start
                    data[itera + 1, 7, 1] = np.array((0.0, 40.301e3))        #Strong kick as it needs to go all the way to Jupiter

                if itera*self.deltat >= 1.96e8 and n <= 1:        #Just a little loop to display the velocity of the Jovian satellite
                    n += 1                                      #once after slingshot to see the speed in relation to the speed of light
                    print("Here is the speed of the Jovian satellite: {0:6.4}".format(self.modulus(data[itera, 7, 1]) / 300000)+"c")






                #Side note, I did this ratio because I actually managed to make a slingshot going at a altitude lower
                #than Jupiter's radius and propel to 33.79 the speed of light. Voilà

    def euler(self, planet, pos, vel):      #Method to calculate the position and Velocity using the Euler-Cromer method

        a = self.wholedata[0, planet, 2]    #Import the acceleration of the current iteration
        vel = vel + a*self.deltat
        pos = pos + vel*self.deltat
        return pos, vel         #Return the newly calculated position and velocity to be appended in the dataset

    def accelcalc(self, po1, po2, m):       #Method for calculating the acceleration of a body using both position vector
                                            #of the planet of interest and the other body
        rvec = po1 - po2    #Create the actual vector (distance) between the two bodies by substracting their position vector
        accel = (m/((self.modulus(rvec))**3))*rvec  #Uses the equation in the index of the report, also calls the modulus method
        return accel        #Return the calculated acceleration for furture purpose


    def modulus(self, rvec):        #Calculate the modulus of the given vector

        modul = (rvec[0]**2 + rvec[1]**2)**0.5
        return modul

    def closeapproachjupiter(self):     #Method for calculating and printing the closest distance achieved by satellite to Jupiter

        for i in range(self.itera):     #For loop which will go through the data to calculate the actual vector between
            vecpo = self.wholedata[i, 5, 0] - self.wholedata[i, 7, 0]#the satellite and Jupiter
            modu = self.modulus(vecpo)      #Calculate the modulus of it to get the actual distance
            if modu < 10e8:     #If this distance is lower than 1,000,000 km then:
                a = i * self.deltat     #Mulitply i by the timestep to get the acual time at which it occured
                approach = a / 60.0 / 60.0 / 24.0 / 365.0       #Divide it by the numer of seconds in a year on Earth to get it in an understandable way
                print("The closest distance to Jupiter: {0:10.2}".format(modu) + " m, and at {0:5.3}".format(approach)+" years")

    def closeapproachmars(self):        #Similar mether than above but for Mars

        for i in range(self.itera):
            vecpo = self.wholedata[i, 4, 0] - self.wholedata[i, 6, 0]   #Instead of Jupiter (index 5) it is Mars (index 4)
            modu = self.modulus(vecpo)
            if modu < 10e6:     #Lower distance as a precise approach was desired
                a = i * self.deltat
                approach = a/60.0/60.0/24.0/365.0
                print("The closest distance to Mars: {0:9.2}".format(modu) + " m, and at {0:5.3}".format(approach)+" years")

    def kinetic(self):      #Method for calculating the total energy of system
        data = self.wholedata       #Again, self.wholedata will be used extensively so the short "data" is used for convenience
        ktime = int(self.itera/100.0)       #create a integer which represents the number of times the energy was calculated
        self.tkedata = np.zeros(ktime)      #create a numpy array of this size to store the total energy
        self.n = 0      #Creates a counter which will be used for plotting as the x axis values
        for i in range(self.itera):     #For loop to go through every iterations
            mod = i % 100     #Creates a modulo of i
            if mod == 0:      #Check when the modulo is equal to zero, meaning 100 iterations has passed
                self.n += 1     #Increase the counter by 1
                tke = 0        #Reinitialise the tke (sum holder for the total energy
                for j in range(self.numlines):      #For loop which go through all planets

                    ke = 0.5 * self.masslist[j] * (self.modulus(data[i, j, 2]))**2      #Calculates kinetic energy using the
                                                                                    #Equation in the report
                    tke += ke               #Adds the kinetic energy to the total energy
                    for k in range(self.numlines):      #for loop to go through all the Other planets
                        rvec = data[i, j, 0] - data[i, k, 0]    #Calculates the vector position between the two bodies

                        if j != k:      #When the two planets are different
                            pge = -self.G * ((self.masslist[k] * self.masslist[j]) / self.modulus(rvec))#Calculate the
                            #potential gravitational energy using the equation in the report
                            tke += pge  #Adds it to the total energy

                self.tkedata[self.n - 1] += tke     #Append the total energy in the data holder to later plot it.

    def orbitalperiod(self):    #Method for calculating the orbital period for the planets.

        data = self.wholedata
        for i in range(1, self.numlines-2): #for loop for going through all bodies
                # #Limit the application of this loop to only the planets (thus excluding the Sun and the satellites
                orbit = math.sqrt((data[0, i, 0, 0]**3 * 4*math.pi**2)/(self.G*(self.masslist[i]+self.masslist[0])))
                #Calculate the time of the orbit using derived equation of Third law Kepler (in the report)
                orbit = float(orbit/60.0/60.0/24.0/365.0)#Adapt it into Earth years
                print("The orbital period of " + str(self.bodyname[i]) + " is : {0:8.4}".format(orbit) + " Earth years")
                #Print it the its name.

    def animate(self, i):       #Short method, but crucial for the Funcanimation library
        for j in range(self.numlines):      #For loop which put the centre of the patches using the position coordinates
            self.patchli[j].center = (self.wholedata[i, j, 0, 0], self.wholedata[i, j, 0, 1])

        return self.patchli  # No coma as it is iterable

    def display(self):      #Method used for displaying the dataset
        fig = plt.figure()
        ax1 = fig.add_subplot(121)      #Divide the plot into two subplot. ax1 will be the animated solar system
        ax2 = fig.add_subplot(122)      #and ax2 will be the total energy plot

        self.patchli = []       #Creates the list which will hold the list of patches

        self.patchli.append(plt.Circle((0, 0), radius=4 * 10 ** 10, color='r', animated=True))  #Append the Sun at 0,0 coordinates (centre)
        for j in range(1, self.numlines):   #For loop for appending the rest of the bodies
            if j <= 5:  #This is for the planets
                self.patchli.append(plt.Circle((self.wholedata[0, j, 0, 0], 0), radius=2e10, color='g', animated=True))
            if j >= 6:  #This is for the satellites, this permits easier visual recognition with different size and colour
                self.patchli.append(plt.Circle((self.wholedata[0, j, 0, 0], 0), radius=1e10, color='b', animated=True))

        for i in range(self.numlines):  #Add the patches to the axe
            ax1.add_patch(self.patchli[i])

        ax1.axis('scaled')
        ax1.set_xlim(-8e11, 8e11), ax1.set_ylim(-8e11, 8e11)       #Set the limits for the x and y axis
        ax1.patch.set_facecolor('black')        #Change the in black because it is space
        ax1.set_title("Solar system animation") #Set a title

        anim = FuncAnimation(fig, self.animate, frames=self.itera, interval=0.1, blit=True) #Function to animate it

        planetplot = np.zeros((self.numlines, 2, self.itera))       #Creates a empty array for plotting the path of planets
#the first dimension will be the planets. For each planet there will be two arrays which contains the position in at every iterations
        for n in range(self.itera):     #For loop to go through iterations
            for j in range(self.numlines):  #For loop to go through planets
                planetplot[j, 0, n] = self.wholedata[n, j, 0, 0]    #Appends the x value of vector position of planet "j" at time "n"
                planetplot[j, 1, n] = self.wholedata[n, j, 0, 1]    #Appends the y value of vector position of planet "j" at time "n"
        for i in range(self.numlines):
            ax1.plot(planetplot[i, 0], planetplot[i, 1])    #plot this data onto the plot on which the animation happens

        ax2.set_xlim(0, int(self.itera / 100.0)), ax2.set_ylim(-6.8e35, -6.6e35)    #Sets the limits for the energy plot
        ax2.set_title("Total energy plot")  #Put a title
        ax2.plot(range(self.n), self.tkedata)   #Plot the actual data

        plt.show()      #Shows the magic of programming to the world
