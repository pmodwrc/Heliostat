from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

import time

import numpy as np
import can
import can.interfaces.kvaser


class IO(QObject):
    """
    Class representing Input/Output operations for controlling motors and reading encoder values via a CAN Bus.

    Attributes:
        fourq (object): An object representing the class that handles the communication with the 4Q Sensor
        A (int): Last known Azimuth position.
        E (int): Last known elevation position.
        data (dict): Dictionary for position data and data from the 4Q sensor
        FQ_out (dict): Dictionary for storing raw output data from the 4Q sensor.
        motor_ID_A (int): ID for Azimuth motor.
        motor_ID_E (int): ID for Elevation motor.
        motor_reply_ID_E (int): Reply ID for Elevation motor.
        motor_reply_ID_A (int): Reply ID for Azimuth motor.
        enc_ID_E (int): ID for Elevation encoder.
        enc_ID_A (int): ID for Azimuth encoder.
        enc_SID_E (int): Elevation encoder SID.
        enc_SID_A (int): Azimuth encoder SID.
        N_no_answer (dict): Dictionary to track no answer counts.
        good_to_send (dict): Dictionary to track send statuses.
        status_reply (dict): Dictionary to track reply statuses.
        M_status (dict): Dictionary to track motor statuses.
        M_rep_C (dict): Dictionary to track motor reply C.
        M_rep_M (dict): Dictionary to track motor reply M.
        Enc_offset_A (int): Offset for Azimuth encoder.
        Enc_offset_E (int): Offset for Elevation encoder.
        scaler_enc (float): Encoder scaler.
        scaler_step (int): Step scaler.
        drehmatrix (list): Rotation matrix.
        bus (can.interface.Bus): CAN bus interface.
    """
    def __init__(self,fourq):         #fourq einfügen
        """
        Initialize the IO class.

        Args:
            fourq (object): An object representing the class that handles the communication with the 4Q Sensor
        """ 
        self.fourq = fourq
        self.A= 180            # last known Azimuth position
        self.E= 63              # last known elevation position
        
        self.data = {}
        self.FQ_out = {}
        self.motor_ID_A = 1
        self.motor_ID_E = 3
        self.motor_reply_ID_E = 4
        self.motor_reply_ID_A = 2
        self.enc_ID_E = 0x185  
        self.enc_ID_A = 0x184
        self.enc_SID_E = 1541
        self.enc_SID_A = 1540          
        
        self.N_no_answer = {'E': 0, 'A': 0}
        self.good_to_send = {'E': 1, 'A': 1}
        self.status_reply = {'E': 0, 'A': 0}
        self.M_status = {'E': 100, 'A': 100}
        self.M_rep_C = {'E': 0, 'A': 0}
        self.M_rep_M = {'E': 0, 'A': 0}
        
        self.Enc_offset_A=141
        
        self.Enc_offset_E=-217
        
        self.soft_end_switch = 0
        
        self.scaler_enc=2**18/-360
        self.scaler_step=1
        
        self.drehmatrix = [[0,-1],[1,0]]
        
        
        self.bus = can.interface.Bus(bustype="kvaser", channel=0, bitrate=500000)
        
        
        self.bus.set_filters([{"can_id": 0, "can_mask": 0x7FF, "extended": False},
                              {"can_id": 2, "can_mask": 0x7FF, "extended": False}])
        self.bus.set_filters([{"can_id": self.enc_ID_A, "can_mask": 0x1FFFFFFF, "extended": True},
                              {"can_id": self.enc_ID_E, "can_mask": 0x1FFFFFFF, "extended": True},
                              {"can_id": self.motor_reply_ID_E, "can_mask": 0x1FFFFFFF, "extended": True},
                              {"can_id": self.motor_reply_ID_A, "can_mask": 0x1FFFFFFF, "extended": True}])
        
        
        # Initialize Encoders and Motor End Switches
        self.initialize_encoder(self.enc_SID_A)
        self.initialize_encoder(self.enc_SID_E)
        
        #Initialize motor end switch
        self.initialize_motor_end_switches()
        
        # Disable Power Off at Zero Velocity for Motor
        self.disable_power_off_at_zero_velocity()
        
    
       
    def send_data(self, M):
        """
        Send data to motors based on input values.

        Args:
            M (dict): A dictionary containing speed values for the motor.

        Returns:
            None
        """
        
        if M['E'] <= 0:
            direction = 1       #rotate right
            speed = int(abs(M['E']*self.scaler_step))
            byte0 = direction.to_bytes(1, "big")
            byte1 = speed.to_bytes(6, "big")
            byte= byte0+byte1
            self.msg_E = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=byte)   
        else:
            direction = 2       #rotate left
            speed = int(abs(M['E']*self.scaler_step))
            byte0 = direction.to_bytes(1, "big")
            byte1 = speed.to_bytes(6, "big")
            byte= byte0+byte1
            self.msg_E = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=byte)
    
            
        
        
        if M['A'] <= 0:
            direction = 1       #rotate right
            speed = int(abs(M['A']*self.scaler_step))
            byte0 = direction.to_bytes(1, "big")
            byte1 = speed.to_bytes(6, "big")
            byte= byte0+byte1
            self.msg_A = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_A,is_fd=False,data=byte)    
        else:
            direction = 2       #rotate left
            speed = int(abs(M['A']*self.scaler_step))
            byte0 = direction.to_bytes(1, "big")
            byte1 = speed.to_bytes(6, "big")
            byte= byte0+byte1
            self.msg_A = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_A,is_fd=False,data=byte)  
        
        

        
        
        if self.good_to_send['E'] == 1:
            #print(self.msg_E)
            self.bus.send(self.msg_E)
            self.good_to_send['E'] = 0
        else:
            print('no reply from motor E')
            self.N_no_answer['E'] += 1
            if self.N_no_answer['E'] > 4:
                self.N_no_answer['E'] = 0
                self.good_to_send['E'] = 1

        if self.good_to_send['A'] == 1:
            #print(self.msg_A)
            self.bus.send(self.msg_A)
            self.good_to_send['A'] = 0
            
        else:
            print('no reply from motor A')
            self.N_no_answer['A'] += 1
            if self.N_no_answer['A'] > 4:
                self.N_no_answer['A'] = 0
                self.good_to_send['A'] = 1
        return self.soft_end_switch
    
    def getdata(self):
        """
        Get data from motors that are already processed.

        Returns:
            dict: A dictionary containing motor data and positional data from the 4Q Sensor.
        
        """
        self.data['A'] = self.A
        self.data['E'] = self.E
        self.data['FQ'] = self.get_4Q_data()  
        return self.data
    
    def get_4Q_data(self):
        """ 
        Gets the data from the 4Q sensor via the Q4_com_module

        Returns:
            numpy.ndarray: Transformed data.
        """
        self.fq = self.fourq.get_data()
        self.FQ_out['X'] = self.fq['X']
        self.FQ_out['Y'] = self.fq['Y']
        self.t1 = np.dot(self.drehmatrix, np.array([self.fq['X'],self.fq['Y']]))
        self.variable_drehmatrix = [[np.cos(-self.A),-np.sin(-self.A)],[-np.sin(-self.A),-np.cos(-self.A)]]
        
        self.t2 = np.dot(self.variable_drehmatrix,self.t1)
        return self.t2
        
        
        
    def request_manual_input(self):
        """
        Request manual input for motors to controll the Heliostat with the switches on the outside of the building.
        
        Returns:
            None
        """
        data = b'\x0f\xff\x00\x00\x00\x00\x00'
        self.msg_E = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(self.msg_E)
        
        data = b'\x0f\xff\x00\x00\x00\x00\x00' 
        self.msg_A = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_A,is_fd=False,data=data)
        self.bus.send(self.msg_A)  
        
   
    def request_motor_data(self):
        """
        Request motor data.
        
        Returns:
            None
        """
        
        
        data = b'\x06\xce\x00\x00\x00\x00\x00'
        self.msg_E = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(self.msg_E)
        
        '''
        data = b'\x06\xce\x00\x00\x00\x00\x00' 
        self.msg_A = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_A,is_fd=False,data=data)
        self.bus.send(self.msg_A)
        '''
        
        data = b'\x06\x0b\x00\x00\x00\x00\x00' 
        self.msg_A = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_A,is_fd=False,data=data)
        self.bus.send(self.msg_A)
        
        
        
        
        
    
    
    def start_notifier(self):
        """
        Start the CAN bus notifier to receive messages asynchronously.
        
        Returns:
            None
        """
        bus = can.interface.Bus(bustype="kvaser", channel=0,bitrate=500000)
        can.Notifier(bus,[self.receive_fcn()])
        
    def receive_fcn(self):
        """
        Function to receive CAN messages and process them.

        This function continuously listens for CAN messages and updates motor data and statuses accordingly.
        
        Returns:
            None

        Raises:
            can.CanError: If an error frame is detected during message reception.
        """
        bus = can.interface.Bus(bustype="kvaser", channel=0,bitrate=500000)
        try:
            for msg in bus:
                if msg.arbitration_id == self.enc_ID_A:
                    relevant_bytes = msg.data[:4]
                    self.A = -(int.from_bytes(relevant_bytes, byteorder='little', signed=True))/self.scaler_enc + self.Enc_offset_A
                    if self.A > 360:
                        self.A = self.A - 360
                elif msg.arbitration_id == self.enc_ID_E:
                    relevant_bytes = msg.data[:4]
                    self.E = -(int.from_bytes(relevant_bytes, byteorder='little', signed=True))/self.scaler_enc + self.Enc_offset_E
                    
                elif msg.arbitration_id == self.motor_reply_ID_E:
                    self.M_status['E'] = int.from_bytes(msg.data[1:2], byteorder='little', signed=True)
                    axis_parameter=int.from_bytes(msg.data[2:3], byteorder='little', signed=True)
                    
                    
                    if axis_parameter == 2 or axis_parameter == 15:
                        self.M_rep_C['E'] = axis_parameter
                        self.M_rep_M['E'] = int.from_bytes(msg.data[6:10], byteorder='little', signed=True)/self.scaler_step
                    elif axis_parameter == 6:
                        self.status_reply['E'] = int.from_bytes(msg.data[6:10], byteorder='little', signed=True)
                    
                    
                    if not (self.M_status['E'] & 100):
                        print("Data Transmission Error Motor E")
                    else:
                        self.good_to_send['E'] = 1
                        
                       
                elif msg.arbitration_id == self.motor_reply_ID_A: 
                    self.M_status['A'] = int.from_bytes(msg.data[1:2], byteorder='little', signed=True)
                    axis_parameter = int.from_bytes(msg.data[2:3], byteorder='little', signed=True)
                    
                    
                    if axis_parameter == 2 or axis_parameter == 15:
                        self.M_rep_C['A'] = axis_parameter
                        self.M_rep_M['A'] = int.from_bytes(msg.data[6:10], byteorder='little', signed=True)/self.scaler_step
                    
                    
                    if axis_parameter == 6:
                        self.status_reply['A'] = int.from_bytes(msg.data[6:10], byteorder='little', signed=True)
                        if self.status_reply['A'] == 0 or self.status_reply['A'] == 1:
                            self.status_soft_end_switch = int.from_bytes(msg.data[6:10], byteorder='little', signed=True)
                    if not (self.M_status['A'] & 100):
                        print("Data Transmission Error Motor E")
                    else:
                        self.good_to_send['A'] = 1
                        
        except can.CanError:
            print("Error frame detected")

        
    
    #Functions to initialize the motors and the encoder
    
    
    def initialize_encoder(self, enc_SID):
        """
        Initialize an encoder.

        Args:
            SID (int): Encoder SID.
        """
        messageout1 = can.Message(arbitration_id=enc_SID, is_extended_id=False, dlc=8)
        messageout1.data = b'\x22\x00\x18\x02\xfe\x00\x00\x00'
        self.bus.send(messageout1)
        
        messageout2 = can.Message(arbitration_id=enc_SID, is_extended_id=False, dlc=8)
        messageout2.data = b'\x22\x00\x18\x05\xe8\x03\x00\x00'
        self.bus.send(messageout2)
        
        messageout3 = can.Message(arbitration_id=0, is_extended_id=False, dlc=8, data=[0] * 8)
        node_id = [(enc_SID - 0x600) & 0xFF]
        node_id = node_id[0]
        data = bytes([1,node_id,0,0,0,0,0,0])
        messageout3.data = data
        self.bus.send(messageout3)
        
    def initialize_motor_end_switches(self):
        """
        Disable motor end switches.

        Args:
            motor_ID (int): Motor ID.
        """
        # Left end switch
        messageout = can.Message(arbitration_id=self.motor_ID_E, is_extended_id=True, dlc=7)
        messageout.data = b'\x05\x0d\x00\x00\x00\x00\x01'
        self.bus.send(messageout)

        # Right end switch
        messageout.data = b'\x05\x0c\x00\x00\x00\x00\x01'
        self.bus.send(messageout)
        
        # Left end switch
        messageout = can.Message(arbitration_id=self.motor_ID_A, is_extended_id=True, dlc=7)
        messageout.data = b'\x05\x0d\x00\x00\x00\x00\x00'
        self.bus.send(messageout)

        # Right end switch
        messageout.data = b'\x05\x0c\x00\x00\x00\x00\x00'
        self.bus.send(messageout)
        

    def disable_power_off_at_zero_velocity(self):
        """
        Disable power off at zero velocity for a motor. !!!Attention with this settings!!! Study the manual before changing the power for the motor.

        Args:
            motor_ID (int): Motor ID.
        """
        messageout = can.Message(arbitration_id=self.motor_ID_E, is_extended_id=True, dlc=7)
        messageout.data = b'\x05\x07\x00\x00\x00\x00\x10'
        self.bus.send(messageout)
        
    def stop_encoder(self):
        """
        Sends stop commands to the encoder.

        Returns:
            None
        """
        data = b'\x02\x04\x00\x00\x00\x00\x00\x00'
        messageout1 = can.Message(is_extended_id=False,dlc=8,arbitration_id=0,is_fd=False,data=data)
        self.bus.send(messageout1)
        
        data = b'\x02\x05\x00\x00\x00\x00\x00\x00'
        messageout2 = can.Message(is_extended_id=False,dlc=8,arbitration_id=0,is_fd=False,data=data)
        self.bus.send(messageout2)
    
    
    
    
    def rotate_right(self):
        """
        This turns the motor to the right infinitely long
        """
        while True:
            vel = 5000
            direction = 1       #rotate right
            speed = int(abs(vel))
            byte0 = direction.to_bytes(1, "big")
            byte1 = speed.to_bytes(6, "big")
            byte= byte0+byte1
            print(byte)
            self.msg_E = can.Message(is_extended_id=True,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=byte)
            self.bus.send(self.msg_E)
            time.sleep(1)
            print("send")
    
        
    def getting_params(self):
        print("I am here")
        data = b'\x06\x04\x00\x00\x00\x00\x00'
        messageout1 = can.Message(is_extended_id=False,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(messageout1)
        print("send 1")
        time.sleep(5)
        
        data = b'\x06\x05\x00\x00\x00\x00\x00'
        messageout1 = can.Message(is_extended_id=False,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(messageout1)
        print("send 2")
        time.sleep(5)
        
        data = b'\x06\x06\x00\x00\x00\x00\x00'
        messageout1 = can.Message(is_extended_id=False,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(messageout1)
        print("send 3")
        time.sleep(5)
        
        data = b'\x06\x07\x00\x00\x00\x00\x00'
        messageout1 = can.Message(is_extended_id=False,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(messageout1)
        print("send 4")
        time.sleep(5)
        
        data = b'\x06\x8c\x00\x00\x00\x00\x00'
        messageout1 = can.Message(is_extended_id=False,dlc=7,arbitration_id=self.motor_ID_E,is_fd=False,data=data)
        self.bus.send(messageout1)
        print("send 5")
        time.sleep(5)
        
        
    def delete(self):
        """
        Calls the stop function of the encoder and shuts down the CAN bus interface
        
        Returns:
            None
        """
        M = {'A': 0, 'E': 0}
        self.send_data(M)
        self.stop_encoder()
        self.bus.shutdown()
        del self.bus