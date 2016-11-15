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
        self.initA() # Reinitialize A
        self.setH() # Update h

    def setH(self):
        self.h = self.length / float(self.n)

    def setE(self, E):
        self.E = E

    def setI(self):
        self.I = (self.width*pow(self.depth, 3)/12.0)

    def initX(self):
        self.x = []
        temp = 0
        for i in range(self.n):
            temp += self.h # x_n = x_n-1 + h
            self.x.append(temp)
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
        """y(x) represents the correct solution for step 2 where x = L (end of the beam) and f is constant"""
        # y(x) = (f / (24*E*I)x^2(x^2 - 4*L*x + 6*L^2)
        return (self.fConst() / (24.0*self.E*self.I)) * pow(x, 2) * (pow(x, 2) - 4.0*self.length*x + 6.0*pow(self.length, 2))
          
        
    def Activity1(self):
        """Activity 1 - Solve for each Y"""
        # A*y = b , solve for y
        # bi = ( h^4/(E*I) )*f(xi)
        b = np.zeros(shape=(self.n,1)) # Initialize b as an nx1 matrix with each entry set to zero
        bi = (pow(self.h, 4) / ( self.E * self.I ))*self.fConst() # Calcualte each entry in b (assuming f(x1)=f(x2)=...=f(xn)=fConst())

        for i in range(self.n):
            b[i][0] = bi

        self.y = np.linalg.solve(self.A, b)

    def Activity2(self, printer):
        """Activity 2 - Plot Solution from step 1 agains the correct solution and check the error at the end of the beam"""

        # TODO: Plot solution from step 1 with yi = self.y[i][0] and xi = self.x[i] for i = 0,...n-1

        y_L = self.y_x(self.length) # Calculate the correct solution at the end of the beam
        if (printer == 5):
            print 'y(',self.length,') =\t',y_L
            print 'Calculated y:\t',self.y[self.n-1][0]
            print 'Error = \t', abs(self.y[self.n-1][0] - y_L)
        return (abs(self.y[self.n-1][0] - y_L))

    def Activity3(self):
        print 'n \t\t k \t\t error \n'
        for k in range (1, 12):
            changingN = 10 * (pow(2,k))
            EBB = EulerBernoulliBeam(2.0, 0.3, 0.03, changingN )
            EBB.Activity1()
            
            error = EBB.Activity2(0)
            print changingN, '\t\t', k , '\t\t', error
        
        
# Tests

# Activity 1
print 'Activity 1:'
EBB = EulerBernoulliBeam(2.0, 0.3, 0.03, 10)
EBB.Activity1()
print 'Y_i:\n',EBB.y
print '\n'

# Activity 2
print 'Activity 2:'
EBB.Activity2(5)
print '\n'

# Activity 3
print 'Activity 3:'
EBB.Activity3()