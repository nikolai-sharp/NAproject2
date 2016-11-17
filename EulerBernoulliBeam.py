# Numerical Analysis 1
# Group Assignment #2 - The Euler-Bernoulli Beam
#
# Authors:
#   David Andrews
#   Giovanni Arias
#   Cameron Bramwell
#   William Rooney
#   Nikolai Sharp
#   Kaiqi Zhang

import numpy as np
import matplotlib.pyplot as plt

class EulerBernoulliBeam:
    def __init__(self, length, width, depth, n):
        self.debug = False

        self.length = float(length)
        self.width = float(width)
        self.depth = float(depth)
        self.n = n
        self.E = 1.3E10
        self.g = 9.81

        self.setI()
        self.initA()
        self.setH()
        self.initX()

    def setL(self, l):
        self.length = float(l)
        self.setH() # Update h

    def setW(self, w):
        self.width = float(w)
        self.setI() # Update I

    def setD(self, d):
        self.depth = float(d)
        self.setI() # Update I

    def setN(self, n):
        self.n = n
        self.initA()    # Reinitialize A
        self.setH()     # Update h
        self.initX()    # Reinitialize X

    def setH(self):
        self.h = self.length / float(self.n)

    def setE(self, E):
        self.E = E

    def setI(self):
        self.I = (self.width*pow(self.depth, 3)/12.0)

    def initX(self):
        self.x = []
        temp = 0
        for i in range(self.n+1):
            self.x.append(temp)
            temp += self.h # x_n = x_n-1 + h
        if (self.debug): print 'X: ',self.x


    def initA(self):
        """Initialize Matrix A"""

        # Check min step number
        if (self.n < 5):
            print '[EulerBernoulliBeam:initA] Warning: At least 5 grid steps required to initialize matrix A: steps set to 5'
            self.n = 5
            self.setH() # Update h

        self.A = np.zeros(shape=(self.n,self.n)) # initialize an nxn matrix with each entry set to zero

        # Row 1
        self.A[0][0] = 16.0
        self.A[0][1] = -9.0
        self.A[0][2] = 8.0/3.0
        self.A[0][3] = -0.25

        # Row 2
        self.A[1][0] = -4.0
        self.A[1][1] = 6.0
        self.A[1][2] = -4.0
        self.A[1][3] = 1.0

        # Rows 3,...,n-2 (indexes: 2,...,n-3)
        for i in range(2,self.n-2):
            self.A[i][i-2]  = 1.0
            self.A[i][i-1]  = -4.0
            self.A[i][i]    = 6.0
            self.A[i][i+1]  = -4.0
            self.A[i][i+2]  = 1.0

        # Row n-1
        self.A[self.n-2][self.n-4] = 16.0/17.0
        self.A[self.n-2][self.n-3] = -60.0/17.0
        self.A[self.n-2][self.n-2] = 72.0/17.0
        self.A[self.n-2][self.n-1] = -28.0/17.0

        # Row n
        self.A[self.n-1][self.n-4] = -12.0/17.0
        self.A[self.n-1][self.n-3] = 96.0/17.0
        self.A[self.n-1][self.n-2] = -156.0/17.0
        self.A[self.n-1][self.n-1] = 72.0/17.0

        
        #print A
        # Print Matrix A without column alignments (Compressed view) when self.debug is set to True
        if (self.debug):
            print 'Matrix A:'
            for i in range(self.n):
                print '[ ',
                for j in range (self.n-1):
                    print self.A[i][j],', ',
                print self.A[i][self.n-1],' ]'

    def fConst(self):
        """f(x) represents only the weight of the beam itself (i.e. with no payload)"""
        return -480.0*self.width*self.depth*self.g

    def y_x(self, x):
        """y(x) represents the correct solution where f is constant"""
        # y(x) = (f / (24*E*I)x^2(x^2 - 4*L*x + 6*L^2)
        if x == 0:
            return 0
        else:
            return (self.fConst() / (24.0*self.E*self.I)) * pow(x, 2) * (pow(x, 2) - 4.0*self.length*x + 6.0*pow(self.length, 2))
          
    def plotA2(self):
        plt.plot(self.x, self.yCalculated, '-r', label='Calculated')
        plt.plot(self.x, self.yActual, '-b', label='Actual')
        plt.legend(loc='lower left')
        plt.show()

    def plotA6(self):
        plt.plot(self.x, self.yCalculated, '-r')
        plt.show()

    def endError(self):
        """Calculate error at end of beam - Activity1() & calcYActual() must be executed to check error"""
        return abs(self.yCalculated[-1] - self.yActual[-1])

    def calcYActual(self):
        self.yActual = []
        for i in self.x:
            self.yActual.append(self.y_x(i)) # Calculate the correct solution at each x_i

    def Activity1(self):
        """Activity 1 - Solve for each Y"""
        self.yCalculated = [0.0]

        # A*y = b , solve for y
        # bi = ( h^4/(E*I) )*f(xi)
        b = np.zeros(shape=(self.n,1)) # Initialize b as an nx1 matrix with each entry set to zero
        bi = (pow(self.h, 4) / ( self.E * self.I ))*self.fConst() # Calcualte each entry in b (assuming f(x1)=f(x2)=...=f(xn)=fConst())

        for i in range(self.n):
            b[i][0] = bi

        yTemp = np.linalg.solve(self.A, b)

        for i in yTemp:
            self.yCalculated.append(i[0])

    def Activity2(self):
        """Activity 2 - Plot Solution from step 1 agains the correct solution and check the error at the end of the beam"""
        self.calcYActual()
        print 'End of Beam Error:',self.endError()
        self.plotA2()

    def Activity3(self):
        print 'n \t\t k \t\t error \n'
        for k in range (1, 12):
            changingN = 10 * (pow(2,k))
            self.setN(changingN)
            self.Activity1()
            self.calcYActual()
            
            error = self.endError()
            print changingN, '\t\t', k , '\t\t', error
        
    def Activity6(self):
        """Activity 6 - Solve for each Y with a 70kg diver balancing on the last 20cm of the beam"""
        
        # TODO: set optimal n based on results from Activity 5
        self.setN(1280) # ***Temporary - Referenced from http://mason.gmu.edu/~zzerhoun/Math447.rc2.5

        self.yCalculated = [0.0]

        # A*y = b , solve for y
        # bi = ( h^4/(E*I) )*f(xi)
        b = np.zeros(shape=(self.n,1)) # Initialize b as an nx1 matrix with each entry set to zero
        bi = (pow(self.h, 4) / ( self.E * self.I ))*self.fConst() # Calcualte each entry in b (assuming f(x1)=f(x2)=...=f(xn)=fConst() and x < 1.8)
        bi_diver = (pow(self.h, 4) / ( self.E * self.I ))*(self.fConst()-self.g*70.0/0.2) # Calcualte each entry in b (assuming f(x1)=f(x2)=...=f(xn)=fConst() and 1.8 <= x <= 2.0)

        for i in range(self.n):
            if self.x[i+1] >= 1.8 and self.x[i+1] <= 2.0: # Diver present
                b[i][0] = bi_diver
            else:
                b[i][0] = bi

        yTemp = np.linalg.solve(self.A, b)

        for i in yTemp:
            self.yCalculated.append(i[0])

        print 'Deflection at x = L:',self.yCalculated[-1]

        self.plotA6()
        
# Tests

# Activity 1
print 'Activity 1:'
EBB = EulerBernoulliBeam(2.0, 0.3, 0.03, 10)
EBB.Activity1()
print 'Y_i:'
for i in EBB.yCalculated:
    print '[',i,']'
print '\n'

# Activity 2
print 'Activity 2:'
EBB.Activity2()
print '\n'

# Activity 3
print 'Activity 3:'
EBB.Activity3()

# Activity 6
print 'Activity 6:'
EBB.Activity6()