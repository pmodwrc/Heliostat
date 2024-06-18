import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *


class PI(QObject):
    """A class implementing a Proportional-Integral controller.

    Attributes:
        array (list): A list to store past error values.
        P (float): Proportional gain.
        I (float): Integral gain.
        grad2step (float): Conversion factor.
        E (float): Current error value.
        V (float): Current process variable value.
        rest (float): Rest value from previous calculations.
        sum (float): Summation of past errors.
    """
    def __init__(self, N_array, P, I, conversion):
        """Initialize the PI controller.

        Args:
            N_array (int): Length of the array to store past error values.
            P (float): Proportional gain.
            I (float): Integral gain.
            conversion (float): Conversion factor.
        """

        self.array = [0] * N_array
        self.P = P
        self.I = I
        self.grad2step= conversion
        self.E = 0
        self.V = 0
        self.rest = 0
        self.sum = 0
        
    def feed(self, E,V):
        """
        Feed the controller with the current error and process variable values.

        Args:
            E (float): Current error value.
            V (float): Current process variable value.
        """
        self.V = V
        self.E = E
        self.array = np.roll(self.array, 1).tolist()
        self.array[0] = self.E


        
    
    def get_speed(self):
        """
        Calculate and return the control output.
        
        Returns:
            None
        """
        self.sum = self.sum +self.E
        self.M = (self.V+np.mean(self.array)*self.P+self.sum*self.I)*self.grad2step + self.rest
        self.rest = self.M-round(self.M)
        return self.M
        
        
    