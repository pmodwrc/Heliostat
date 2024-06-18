from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from datetime import datetime
import os
import numpy as np

import serial


class Q4(QObject):
    """
    A class representing a control interface for a device connected to COM1.

    Attributes:
        Port_4Q (serial.Serial): Serial port connection settings.
        posA (float): Position A value.
        posB (float): Position B value.
        pos (dict): Dictionary to store position data.
        scaler (float): Scaling factor.
        log_directory (str): Directory to store the log file of 4Q sensor data.
        filename (str): Name of the log file.
    """
    def __init__(self):
        """
        Initializers the comunication with the 4Q sensor and makes a logfile for 4Q data.
        
        Returns:
            None
        """
        self.Port_4Q = serial.Serial(port="COM1", baudrate=19200, timeout=1, stopbits=serial.STOPBITS_ONE,parity= serial.PARITY_NONE,bytesize=serial.EIGHTBITS )
        self.Port_4Q.set_buffer_size(rx_size=10000)
        self.posA = 0
        self.posB = 0
        self.pos = {}
        self.scaler = 2.6*0.0151
        self.log_directory = "log"
        now = datetime.now()
        self.filename = os.path.join(self.log_directory, now.strftime("%Y-%m-%d_%H-%M") + "_Heliostat_pointing.dat")
        with open(self.filename, "w"):
            pass
    
    def go(self):
        """
        Asks for data from the device and transforms it to usable data for the GUI.
    
        Returns:
            None
        """
        self.Port_4Q.write(b's')
        now = datetime.now()
        while True:
            str_4Q = self.Port_4Q.readline().decode().strip()
            if str_4Q:
                break
        try:
            ind1 = str_4Q.find('*')
            values_str = str_4Q[ind1 + 1 : -1]
            out_values_4Q = [int(value) for value in values_str.split(',')]
        except ValueError:
            out_values_4Q = [float('nan')] * 6 
        with open(self.filename, "a") as file:
            now = datetime.now()
            file.write(now.strftime("%H:%M:%S "))
            self.fourQ_out = ' '.join(str(x) for x in out_values_4Q)
            file.write(self.fourQ_out + "\n")
        self.fourQ_out = [int(x) for x in out_values_4Q]
        self.posA = ((self.fourQ_out[0]+self.fourQ_out[2])-(self.fourQ_out[1]+self.fourQ_out[3]))/np.sum(np.array(self.fourQ_out))*self.scaler
        self.posB = ((self.fourQ_out[0]+self.fourQ_out[1])-(self.fourQ_out[2]+self.fourQ_out[3]))/np.sum(np.array(self.fourQ_out))*self.scaler
        
        if np.sum(out_values_4Q) < 2000:
            self.posA = 0
            self.posB = 0   

    def get_data(self):
        """
        Retrieve the position data.

        Returns:
            dict: A dictionary containing the positional data from the 4Q sensor.
            
        """
        
        self.pos['X'] = self.posA
        self.pos['Y'] = self.posB
        return self.pos
    
    def delete(self):
        """ 
        Closes the serial connection to the 4Q sensor
        
        Returns:
            None
        """
        self.Port_4Q.close()
        
        
        
        