""" This file contains the GUI as well as the main calculations for controlling the Heliostat.
Autor: L. Bertoli
Date: 21.05.2024"""

import sys
import os
import webbrowser
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.uic import loadUi
import pyqtgraph as pg
from datetime import datetime, timedelta
from pvlib.solarposition import get_solarposition
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

from IO_module import IO
from PI_module import PI
from PRM_module import PRM
from Q4_com_module import Q4


def resource_path(relative_path):
    """
    Get the absolute path to a resource, considering both development and PyInstaller packaging.

    Args:
        relative_path (str): The relative path to the resource.

    Returns:
        str: The absolute path to the resource.
    """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS2
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

class WorkerReceiving(QObject):
    """
    A class representing a worker for receiving data.

    Attributes:
        io (IO): The IO module to handle CAN Bus communication
    """
    def __init__(self, io):
        """
        Initialize the WorkerReceiving.

        Args:
            io (IO): An instance of the IO class for input/output operations.
        """
        super().__init__()
        self.io = io
    def run(self):
        """
        Run the worker loop for receiving data.

        This method continuously calls the start_notifier method of the IO object to receive data.
        """
        while True:
            self.io.start_notifier()
            
class WorkerLogger(QObject):
    """
    A class representing a worker for logging data.

    Attributes:
        log_file_path (str): Path to the log file.
    """
    def __init__(self, log_file_path):
        """
        Initialize the WorkerLogger.

        Args:
            log_file_path (str): Path to the log file.
        """
        self.log_file_path = log_file_path
        if os.path.exists(self.log_file_path):
            pass
        else:
            with open(self.log_file_path, "w") as file:
                file.write('time state POS_A POS_E FQ_X FQ_Y V_A V_E E_A E_E\n')
    def log(self, time, state, pos, M, E):
        """
        Log data to the log file.

        Args:
            time (float): Time value.
            state (str): State information.
            pos (dict): Dictionary containing position information.
            M (dict): Dictionary containing motor speed values.
            E (dict): Dictionary containing error values.
        """
        with open(self.log_file_path, "a") as file:
            file.write(f"{round(time,4)} {state} {pos['A']} {pos['E']} {pos['FQ'][0]} {pos['FQ'][1]} {M['A']} {M['E']} {E['DA']} {E['DE']}\n")   
        
        
            
                   
class RealTimeData(QWidget):
        
    """
    A class representing a widget for displaying real-time data.

    Attributes:
        named_widget (QWidget): The widget containing real-time data.
        position_label_a (QLabel): Label for position data of axis A.
        position_label_e (QLabel): Label for position data of axis E.
        speed_label_a (QLabel): Label for speed data of axis A.
        speed_label_e (QLabel): Label for speed data of axis E.
        mode_label (QLabel): Label for displaying current mode.

    """
    
    def __init__(self):
        """
        Initialize the RealTimeData widget.
        """
        super(RealTimeData, self).__init__()
        loadUi("ui_files/RealTimeData.ui", self)
        self.named_widget = self.findChild(QWidget, "real_time_data")
        self.position_label_a = self.findChild(QLabel, "position_label_a")
        self.position_label_e = self.findChild(QLabel, "position_label_e")
        self.speed_label_a = self.findChild(QLabel, "speed_label_a")
        self.speed_label_e = self.findChild(QLabel, "speed_label_e")
        self.mode_label = self.findChild(QLabel, "mode_label")
    def update_variables(self,M,pos):
        """
        Update the displayed variables based on the provided data.

        Args:
            M (dict): Dictionary containing speed data for axes A and E.
            pos (dict): Dictionary containing position data for axes A and E,
                        and active sense data for axes F and Q.
        """
        self.position_label_a.setText(str(round(pos['A'],2)))
        self.position_label_e.setText(str(round(pos['E'],2)))
        self.speed_label_a.setText(str(round(M['A'], 2)))
        self.speed_label_e.setText(str(round(M['E'],2)))
        self.active_sense_label_a.setText(str(round(pos['FQ'][0])))
        self.active_sense_label_e.setText(str(round(pos['FQ'][1])))
    def update_mode(self, mode):
        """
        Update the mode label with the provided mode.

        Args:
            mode (str): The mode to be displayed.
        """
        self.mode_label.setText(mode)
    
    
class Diagnose(QWidget):
    """
    A class representing a widget for diagnostic motor parameters.

    Attributes:
        graphWidget_E (PlotWidget): Widget for displaying primary motor speed graph.
        graphWidget_A (PlotWidget): Widget for displaying secondary motor speed graph.
        data_line_E (PlotDataItem): Data line for primary motor speed graph.
        data_line_A (PlotDataItem): Data line for secondary motor speed graph.
        M_E (list): List to store primary motor speed data.
        M_A (list): List to store secondary motor speed data.
    
    """
    
    def __init__(self):
        """
        
        Initialize the Diagnose widget.
        
        """
        super(Diagnose, self).__init__()
        loadUi("ui_files/Diagnose.ui", self)
        self.findChild(QWidget, "diagnose_parameters")
                
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.motor_primary.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_E = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">Motor speed primary [µsteps/s]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.motor_secondary.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_A = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">Motor speed secondary [µsteps/s]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        self.M_E = []
        self.M_A = []
    def update_diagnose_plot(self, t, M):
        """
        Update the motor speed graphs with new data.

        Args:
            t (list): List of time values.
            M (dict): Dictionary containing motor speed data for axes A and E.
        """
        self.M_E.append(M['E'])
        if len(self.M_E) == 100:
            self.M_E = self.M_E[1:]
        self.data_line_E.setData(t,self.M_E)
        
        self.M_A.append(M['A'])
        if len(self.M_A) == 100:
            self.M_A = self.M_A[1:]
        self.data_line_A.setData(t,self.M_A)
        
class FourQ(QWidget):
    """
    A class representing a widget for displaying 4Q position parameters.

    Attributes:
        graphWidget_primary (PlotWidget): Widget for displaying primary axis position graph.
        graphWidget_secondary (PlotWidget): Widget for displaying secondary axis position graph.
        data_line_primary (PlotDataItem): Data line for primary axis position graph.
        data_line_secondary (PlotDataItem): Data line for secondary axis position graph.
        x_axis (list): List to store primary axis position data.
        y_axis (list): List to store secondary axis position data.
    """
   
    def __init__(self):
        """
        
        Initialize the FourQ widget.
        
        """
        super(FourQ, self).__init__()
        loadUi("ui_files/4QParameters.ui", self)
        self.findChild(QWidget, "fourq_parameters")
                
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.primary_axis.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_primary = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">4Q Position primary axis [deg]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.secondary_axis.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_secondary = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">4Q Position secondary axis [deg]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
    
        self.x_axis = []
        self.y_axis = []    
    def update_fourq_plot(self, t, pos):
        """
        Update the position graphs with new data.

        Args:
            t (list): List of time values.
            pos (dict): Dictionary containing 4Q position data for primary and secondary axes.
        """
        self.x_axis.append(pos['FQ'][0])
        if len(self.x_axis) == 100:
            self.x_axis= self.x_axis[1:]
        self.data_line_primary.setData(t,self.x_axis)
        
        self.y_axis.append(pos['FQ'][1])
        if len(self.y_axis) == 100:
            self.y_axis= self.y_axis[1:]
        self.data_line_secondary.setData(t,self.y_axis)
        
        
class MotionParams(QWidget):
    """
    A class representing a widget for displaying motion parameters.

    Attributes:
        graphWidget_pos_primary (PlotWidget): Widget for displaying primary axis position graph.
        graphWidget_error_primary (PlotWidget): Widget for displaying primary axis error graph.
        graphWidget_pos_secondary (PlotWidget): Widget for displaying secondary axis position graph.
        graphWidget_error_secondary (PlotWidget): Widget for displaying secondary axis error graph.
        data_line_pos_primary (PlotDataItem): Data line for primary axis position graph.
        data_line_error_primary (PlotDataItem): Data line for primary axis error graph.
        data_line_pos_secondary (PlotDataItem): Data line for secondary axis position graph.
        data_line_error_secondary (PlotDataItem): Data line for secondary axis error graph.
        pos_primary (list): List to store primary axis position data.
        error_primary (list): List to store primary axis error data.
        pos_secondary (list): List to store secondary axis position data.
        error_secondary (list): List to store secondary axis error data.

    """
    def __init__(self):
        """
        
        Initialize the MotionParams widget.
        
        """
        super(MotionParams, self).__init__()
        loadUi("ui_files/MotionParameter.ui", self)
        self.findChild(QWidget, "motion_parameters")
                
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.pos_primary.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_pos_primary = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">Primary position [deg]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.error_pri.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_error_pri = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">Error passive PRI [deg]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.pos_secondary.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_pos_secondary = self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">Secondary position [deg]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        #Initiating graph         
        self.graphWidget = pg.PlotWidget()
        self.error_sec.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_error_sec= self.graphWidget.plot([], [], pen=pen,symbol="+",symbolSize=10,symbolBrush="b")
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">Error passive SEC [deg]</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">Time [UTC + 1]</span>")
        self.graphWidget.showGrid(x=True, y=True)
        
        self.pos_primary = []
        self.error_pri = []
        self.pos_secondary = []
        self.error_sec = []
        
    def update_motion_plot(self, t, pos, error):
        """
        Update the motion parameters graphs with new data.

        Args:
            t (list): List of time values.
            pos (dict): Dictionary containing position data for primary and secondary axes.
            error (dict): Dictionary containing error data for primary and secondary axes.
        """
        self.pos_primary.append(pos['A'])
        if len(self.pos_primary) == 100:
            self.pos_primary = self.pos_primary[1:]
        self.data_line_pos_primary.setData(t,self.pos_primary)

        self.error_pri.append(error['DA'])
        if len(self.error_pri) == 100:
            self.error_pri = self.error_pri[1:]
        self.data_line_error_pri.setData(t,self.error_pri)
        
        self.pos_secondary.append(pos['E'])
        if len(self.pos_secondary) == 100:
            self.pos_secondary = self.pos_secondary[1:]
        self.data_line_pos_secondary.setData(t,self.pos_secondary)
        
        self.error_sec.append(error['DE'])
        if len(self.error_sec) == 100:
            self.error_sec = self.error_sec[1:]
        self.data_line_error_sec.setData(t,self.error_sec)
        
class SpatialView(QWidget):
    """
    A class representing a spatial view widget.

    Attributes:
        graphWidget (PlotWidget): A PlotWidget for 2D plots.
        fig (Figure): A Figure object for 3D plots.
        canvas (FigureCanvas): A FigureCanvas for displaying the 3D plot.
        axes (Axes3D): Axes for 3D plot.
        prm (PRM): An instance of PRM class.
        a_ref (float): Reference angle in degrees
        home_pos (dict): Dictionary containing home position coordinates.
        ref_circle (numpy.ndarray): Array containing points of reference circle.
        m2 (numpy.ndarray): Array containing points of mirror 2 projected onto the wall of the observatory
        m1 (numpy.ndarray): Array containing points of mirror 1 turned and projected onto the wall of the observatory
    """

    def __init__(self):
        """
        Initialize the SpatialView widget.

        Loads the UI file, initializes 2D and 3D plots, and sets up initial parameters.
        """
        super(SpatialView, self).__init__()
        loadUi("ui_files/SpatialView.ui", self)
        self.findChild(QWidget, "spatial_view")
    
        
        #Initializing 2D plot
        self.graphWidget = pg.PlotWidget()
        self.two_dim.addWidget(self.graphWidget)
       
        pen = pg.mkPen(color=(0, 255, 0))
        self.data_line_ref_circle= self.graphWidget.plot([], [], pen=pen)
        pen = pg.mkPen(color=(0, 0, 255))
        self.data_line_m2_circle= self.graphWidget.plot([], [], pen=pen)
        pen = pg.mkPen(color=(255, 0, 0))
        self.data_line_m1_circle= self.graphWidget.plot([], [], pen=pen)
        
        self.graphWidget.setBackground("w")
        self.graphWidget.setLabel("left", "<span style=\"color:black;font-size:11px\">z - coordinate tangens</span>")
        self.graphWidget.setLabel("bottom", "<span style=\"color:black;font-size:11px\">y - coordinate tangens</span>")
        self.graphWidget.showGrid(x=True, y=True)
        self.graphWidget.setAspectLocked(lock=True)
        
        #Initializing 3D plot
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(111, projection='3d')
        self.three_dim.addWidget(self.canvas)
        self.axes.set_xlim(8000,9000)
        
        
        self.prm = PRM()
        self.a_ref = 2.5
        
        
        

        self.axes.scatter(self.prm.P1c[0],self.prm.P1c[1],self.prm.P1c[2],color='b', marker='o')
        self.axes.scatter(self.prm.P2c[0],self.prm.P2c[1],self.prm.P2c[2], color='b', marker='o')
        
        self.home_pos = {'E':63, 'A':180}
        self.prm.update_position(self.home_pos)
        
        #2D Plots
        self.ref_circle = self.prm.circle_points([0,0,0],[0,0,1], np.tan(np.deg2rad(self.a_ref)), 100)
        x_coords = self.ref_circle[:, 0][1:]
        y_coords = self.ref_circle[:, 1][1:]
        self.data_line_ref_circle.setData(x_coords, y_coords)

        
        self.m2 = self.prm.project_m2()
        x_coords = -(self.m2[:, 0])
        y_coords = self.m2[:, 2]
        self.data_line_m2_circle.setData(x_coords, y_coords)
        
        self.m1 = self.prm.project_m1(self.home_pos)
        x_coords = -(self.m1[:, 0])
        y_coords = self.m1[:, 2]
        self.data_line_m1_circle.setData(x_coords, y_coords)
    
    def plot_vector(self,p1,p2):
        """
        Plot a vector in the 3D plot.

        Args:
            p1 (numpy.ndarray): Starting point of the vector.
            p2 (numpy.ndarray): Ending point of the vector.
        """
        p2 = p2-p1
        self.axes.quiver(p1[0], p1[1], p1[2],p2[0], p2[1], p2[2], normalize=False)
    def plot_line(self, p1, p2, color):
        """
        Plot a line in the 3D plot.

        Args:
            p1 (numpy.ndarray): Starting point of the line.
            p2 (numpy.ndarray): Ending point of the line.
            color (str): Color of the line.
        """
        point1 = p1
        point2 = p2
        self.axes.plot([point1[0], point2[0]], [point1[1], point2[1]], [point1[2], point2[2]], color=color, label='Line')
    def plot_mirror(self, P, n, r, nn,color):
        """
        Plot a mirror in the 3D plot.

        Args:
            P (numpy.ndarray): Center point of the mirror.
            n (numpy.ndarray): Normal vector of the mirror.
            r (float): Radius of the mirror.
            nn (int): Number of points to generate on the mirror.
            color (str): Color of the mirror.
        """
        Cp = self.prm.circle_points(P,n, r, nn)
        x = np.array(Cp[:, 0])[1:]
        y = np.array(Cp[:, 1])[1:]
        z = np.array(Cp[:, 2])[1:]
        self.axes.plot(x,y,z, color = color)
        
        
    def set_axes_equal(self,ax):
        """
        Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc.

        Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
        """

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5*max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius]) 
          
    def update_view_plot(self,pos):
        """
        Update the view plot based on the new position.

        Args:
            pos (dict): Dictionary containing the new position coordinates.

        This method updates the position of the mirrors and other elements in the view plot based on the new position.
        It clears the 3D plot, plots mirrors, normal vectors, vectors from the sun, and lines from M1 to M2 and to the house.
        It also sets the axes equal for better visualization.

        Additionally, it updates the 2D plot by projecting M1.

        """
        self.prm.update_position(pos)
        
        #Updating 3D plot
        self.axes.cla()
        #Plotting the Mirrors
        self.plot_mirror(self.prm.P2c,self.prm.n_p2,self.prm.m2r, self.prm.nn, 'b')
        self.plot_mirror(self.prm.P1c,self.prm.n_p1,self.prm.m2r, self.prm.nn, 'b')
        #Plotting normalvector of mirror 2
        self.plot_line(self.prm.P2c, self.prm.P2c+self.prm.n_p2*1000, 'c')
        #Plotting vector from sun to mirror 1
        self.plot_line(self.prm.P1c, self.prm.P1c+self.prm.V_sun*1000, 'r')
        #Plotting normalvector of mirror 1
        self.plot_line(self.prm.P1c, self.prm.P1c+self.prm.n_p1*1000, 'c')
        #Plotting lines from M1 to M2 and to the house
        self.plot_line(self.prm.P1c,self.prm.P2c,'r')
        self.plot_line(self.prm.PV,self.prm.P2c, 'r')
        self.set_axes_equal(self.axes)
        
        #Updating 2D plot
        self.m1 = self.prm.project_m1(self.home_pos)
        x_coords = -(self.m1[:, 0])
        y_coords = self.m1[:, 2]
        self.data_line_m1_circle.setData(x_coords, y_coords)
        
        
        
                  
         
class HeliostatUI(QMainWindow):
    """
    This class represents the main user interface for controlling a heliostat system.
    It provides functionalities for controlling the heliostat's motion, tracking celestial objects,
    displaying real-time data plots, and accessing diagnostic information.

    Attributes:
        latitude (float): The latitude of the heliostat's location.
        longitude (float): The longitude of the heliostat's location.
        altitude (float): The altitude of the heliostat's location.
        target_ra (float): The right ascension (RA) of the target celestial object.
        target_dec (float): The declination (DEC) of the target celestial object.
        trafoA (float): Transformation factor for azimuth axis.
        trafoE (float): Transformation factor for elevation axis.
        home_pos (dict): Dictionary containing the home position angles for both axes.
        passive_speed (dict): Dictionary containing the passive tracking speeds for both axes.
        M (dict): Dictionary containing the motion parameters for both axes.
        low_limit_A (float): Lower limit for azimuth axis.
        high_limit_A (float): Upper limit for azimuth axis.
        low_limit_E (float): Lower limit for elevation axis.
        high_limit_E (float): Upper limit for elevation axis.
        target_offset (dict): Dictionary containing the offset for the target position for both axes.
        cycle_time (int): Time interval for updating the state in milliseconds.
        M_max (int): Maximum speed for motion parameters.
        M_slow_fixed_target (float): Speed for fixed target tracking.
        M_slow_rel (float): Relative speed for passive tracking.
        M_slow_rel_Az (float): Relative speed for passive tracking in azimuth axis.
        M_fine_rel (dict): Dictionary containing relative speed for fine tracking for both axes.
        M_fine_I (dict): Dictionary containing I value for passive tracking for both axes.
        M_fine_active (dict): Dictionary containing active fine tracking speeds for both axes.
        M_fine_active_I (dict): Dictionary containing I value for active tracking for both axes.
        active_frame (int): Maximum allowed error for active tracking.
        update_period (int): Time interval for updating plots in seconds.
        t0 (float): Initial time for calculating time component.
        timer_count (int): Counter for timer ticks.
        t_vector (list): List containing time values for plotting.
        grad2step (float): Conversion factor from degrees to motor steps.
        E (dict): Dictionary containing error values for both axes.
        diagnose (Diagnose): Instance of Diagnose class for diagnostics plotting.
        motion_params (MotionParams): Instance of MotionParams class for motion parameters plotting.
        q_params (FourQ): Instance of FourQ class for 4Q sensor data plotting.
        spatial_view (SpatialView): Instance of SpatialView class for spatial view plotting.
        real_time_data_window (RealTimeData): Instance of RealTimeData class for real-time data plotting.
        passivPI_E (PI): Instance of PI class for passive tracking in elevation axis.
        passivPI_A (PI): Instance of PI class for passive tracking in azimuth axis.
        activePI_E (PI): Instance of PI class for active tracking in elevation axis.
        activePI_A (PI): Instance of PI class for active tracking in azimuth axis.
        fourq (Q4): Instance of Q4 class for handling 4Q sensor data.
        io (IO): Instance of IO class for input/output operations.
        current_date (str): Current date in "YYYY-MM-DD" format.
        log_file_path (str): Path for log file.
        logger (WorkerLogger): Instance of WorkerLogger class for logging data.
        thread_receiving (QThread): Thread for receiving data from CAN Bus.
        worker_receiving (WorkerReceiving): Worker for receiving data from CAN Bus.
        timer (QTimer): Timer for updating state.
        timer_fourq (QTimer): Timer for updating 4Q data.
    """
    
    def __init__(self):
        """
        
        Initializes the HeliostatUI instance and the main window

        """
        
        # Create Window
        super(HeliostatUI, self).__init__()
        loadUi(resource_path("ui_files/Heliostat.ui"), self)
        self.setWindowTitle("GUI Heliostat")
        
       
        self.state = 0
        self.latitude = 46.81
        self.longitude = 9.84
        self.altitude = 1580
        self.target_ra = 0
        self.target_dec = 0
        
        self.trafoA = -5
        self.trafoE = 27
    
        self.pos = {'A': 180, 'E': 63, 'FQ': [0,0]}
        
        self.home_pos = {'E': 63, 'A': 180}
        
        self.passive_speed = {'E': 0 , 'A': 0}
        self.M = {'E': 0, 'A': 0}
        self.low_limit_A = 83
        self.high_limit_A = 218
        self.low_limit_E = 30
        self.high_limit_E = 91
        
        self.target_offset = {'A': 0, 'E': 0}
        self.target = {'A': 180, 'E':63}
        self.cycle_time = 1000

        self.M_max=5000;  # sps
        self.M_slow_fixed_target=0.1
        self.M_slow_rel=0.07;
        self.M_slow_rel_Az=0.1;
        
        self.M_fine_rel = {'A':-0.0040000000000000036 , 'E': 0}
        self.M_fine_I = {'A': 0 , 'E': 0}
        
        #self.M_fine_relA=-0.1%-0.008; # scaler for fine tracking
        #self.M_fine_relE=-0.01%-0.01; # scaler for fine tracking
        #self.M_fine_I_E= -0.000; # I value for passive tracking 
        #self.M_fine_I_A= -0.000; # I value for passive tracking
        
        self.M_fine_active = {'A': 0.04 , 'E': 0.02}
        self.M_fine_active_I = {'A': 0.0044 , 'E': 0.0004}
        
        #self.M_fine_active_E=0.02 
        #self.M_fine_active_A=0.04
        #self.M_fine_active_I_E=0.0004 
        #self.M_fine_active_I_A=0.0004
        
        self.active_frame = 3
        self.update_period = 2
        time_component = datetime.now()-datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
        self.t0 = time_component.total_seconds()/3600
        self.timer_count = 0
        self.t_vector = []
        
        self.grad2step=256/3.6*160
        self.E = {'DA':0,'DE':0}
        
        #Initializing classes for further use
        self.diagnose = Diagnose()
        self.motion_params = MotionParams()
        self.q_params = FourQ()
        self.spatial_view = SpatialView()
        self.real_time_data_window = RealTimeData()
        self.passivPI_E = PI(10, self.M_fine_rel['E'], self.M_fine_I['E'], self.grad2step)
        self.passivPI_A = PI(4, self.M_fine_rel['A'], self.M_fine_I['A'], self.grad2step)
        self.activePI_E = PI(10, self.M_fine_active['E'], self.M_fine_active_I['E'], self.grad2step)
        self.activePI_A = PI(2, self.M_fine_active['A'], self.M_fine_active_I['A'], self.grad2step)
        
        self.fourq = Q4()
        self.io = IO(self.fourq)         
        self.current_date = datetime.now().strftime("%Y-%m-%d")
        self.log_file_path = f"log/heliostat_log_type2_{self.current_date}.dat"
        self.logger = WorkerLogger(self.log_file_path)
        
        
        #Initializing Thread to receive data from CAN Bus
        self.thread_receiving = QThread()
        self.worker_receiving = WorkerReceiving(self.io)
        self.worker_receiving.moveToThread(self.thread_receiving)
        self.thread_receiving.start()
        self.thread_receiving.started.connect(self.worker_receiving.run)
        
        #Initializing timer to update state
        self.timer = QTimer()
        self.timer.setInterval(self.cycle_time)
        self.timer.timeout.connect(self.update_state)
        self.timer.start()
        
        
        #Initializing timer to get 4Q data
        self.timer_fourq= QTimer()
        self.timer_fourq.setInterval(self.cycle_time)
        self.timer_fourq.timeout.connect(lambda: self.fourq.go())
        self.timer_fourq.start()
        
    
        
        
        
        # Connect GUI widgets to backend 
        self.stop_button.clicked.connect(self.toggle_stop)
        self.home_button.clicked.connect(self.toggle_home)
        self.passiv_tracking_button.clicked.connect(self.toggle_passive_tracking)
        self.active_tracking_button.clicked.connect(self.toggle_active_tracking)
        self.go_to_button.clicked.connect(self.toggle_go_to)
        self.passiv_star_tracking_button.clicked.connect(self.toggle_passive_star_tracking)
        self.manual_remote_button.clicked.connect(self.toggle_manual_remote)
        self.camera_button.clicked.connect(self.toggle_camera)
        
        
        self.real_time_data_button.clicked.connect(self.toggle_real_time_data_window)
        self.diagnose_button.clicked.connect(self.toggle_diagnose)
        self.motion_params_button.clicked.connect(self.toggle_motion_params)
        self.q_params_button.clicked.connect(self.toggle_q_params)
        self.view_button.clicked.connect(self.toggle_view)
        
        self.input_pos_a.valueChanged.connect(self.set_target_pos_a)
        self.input_pos_e.valueChanged.connect(self.set_target_pos_e)
        self.input_active_frame.valueChanged.connect(self.set_active_frame)
        self.input_right_ascension.valueChanged.connect(self.set_ra)
        self.input_declination.valueChanged.connect(self.set_dec)
        self.input_speed.valueChanged.connect(self.set_speed)
       
        
        self.set_target_pos_a(float(self.input_pos_a.text().replace(",", ".")))
        self.set_target_pos_e(float(self.input_pos_e.text().replace(",", ".")))
        self.set_active_frame(float(self.input_active_frame.text().replace(",", ".")))
        self.set_ra(float(self.input_right_ascension.text().replace(",", ".")))
        self.set_dec(float(self.input_declination.text().replace(",", ".")))
        self.set_speed(float(self.input_speed.text().replace(",", ".")))
        
        self.stop_button.setEnabled(True)
        self.home_button.setEnabled(True)
        self.passiv_tracking_button.setEnabled(True)
        self.active_tracking_button.setEnabled(True)
        self.go_to_button.setEnabled(True)
        self.passiv_star_tracking_button.setEnabled(True)
        self.manual_remote_button.setEnabled(True)
        self.camera_button.setEnabled(True)
        
        self.real_time_data_button.setEnabled(True)
        self.diagnose_button.setEnabled(True)
        self.motion_params_button.setEnabled(True)
        self.q_params_button.setEnabled(True)
        self.view_button.setEnabled(True)
    
    def toggle_real_time_data_window(self):
        """
        Toggles for the real time data plot.

        Returns:
            None
        """
        if self.real_time_data_window.isVisible():
            self.real_time_data_window.close()
        else:
            self.real_time_data_window.show()
    
    def toggle_diagnose(self):
        """
        Toggles for the diagnose data plot.

        Returns:
            None
        """
        if self.diagnose.isVisible():
            self.diagnose.close()
        else:
            self.diagnose.show()
        
    def toggle_motion_params(self):
        """
        Toggle for the motion parameters plot.

        Returns:
            None
        """
        if self.motion_params.isVisible():
            self.motion_params.close()
        else:
            self.motion_params.show()
    def toggle_q_params(self):
        """
        Toggle for the 4Q sensor data plot.

        Returns:
            None
        """
        if self.q_params.isVisible():
            self.q_params.close()
        else:
            self.q_params.show()
            
    def toggle_view(self):
        """
        Toggle for the 3D/2D view.

        Returns:
            None
        """
        if self.spatial_view.isVisible():
            self.spatial_view.close()
        else:
            self.spatial_view.show()
    
            
    def set_target_pos_a(self, value):
        """
        Set the target position for the primary axis

        Args:
            value (float): The position to set for the primary axis

        Returns:
            None
        """
        self.target['A'] = value
        
    def set_target_pos_e(self, value):
        """
        Set the target position for the secondary axis

        Args:
            value (float): The position to set for the second axis

        Returns:
            None
        """
        self.target['E'] = value
        
    def set_active_frame(self, value):
        """
        Set the maximum aloud error for active tracking

        Args:
            value (float): The maximum error to set for active tracking

        Returns:
            None
        """
        self.active_frame = value
    
    def set_ra(self, value):
        """
        Set the right ascension (RA) target value.

        Args:
            value (float): The value to set as the target RA.
            
        Returns:
            None
        """
        self.target_ra = value
        
    def set_dec(self, value):
        """
        Set the declination (DEC) target value.

        Args:
            value (float): The value to set as the target DEC.
        
        Returns:
            None
        """
        self.target_dec = value
        
    def set_speed(self, value):
        """
        Set the speed for manual remote controll

        Args:
            value (float): The speed to set for the manual remote mode

        Returns:
            None
        """
        self.speed = value
        
    
    
        
    def toggle_stop(self):
        """
        Toggle to change the current mode to stop

        Returns:
            None
        """
        self.mode.setText("Stop")
        self.state = 0
    def toggle_home(self):
        """
        Toggle to change the current mode to home

        Returns:
            None
        """
        self.mode.setText("Home")
        self.state = 1
    def toggle_go_to(self):
        """
        Toggle to change the current mode to going to a certain position

        Returns:
            None
        """
        self.mode.setText("Fixed position")
        self.state = 2
    def toggle_passive_tracking(self):
        """
        Toggle to change the current mode to passive tracking

        Returns:
            None
        """
        self.mode.setText("Passive tracking")
        self.state = 3
    def toggle_active_tracking(self):
        """
        Toggle to change the current mode to active tracking

        Returns:
            None
        """
        self.mode.setText("Active tracking")
        self.state = 4
    def toggle_passive_star_tracking(self):
        """
        Toggle to change the current mode to star tracking

        Returns:
            None
        """
        self.mode.setText("Passive star tracking")
        self.state = 5
    def toggle_manual_remote(self):
        """
        Toggle to change the current mode to manual remote

        Returns:
            None
        """
        self.mode.setText("Manual remote")
        self.state = 6
    
    def toggle_camera(self):
        """
        Toggle to open the site with access to the camera that is installed remotly. 

        Returns:
            None
        """
        url = 'http://srv-wrc-syno2.ad.pmodwrc.ch:5000/index.cgi?launchApp=SYNO.SDS.SurveillanceStation#/signin'
        webbrowser.open(url)
        
    
    
    def update_state(self):
        """
        This function updates the state of the system according to the current state value.
        It performs various actions based on the state, such as stopping motion, going to the home position,
        going to a fixed target, passive sun tracking, active sun tracking, star tracking, manual remote control,
        and logging data. It also updates real-time data plots, sends speed data to the motor
        and performs limit checks to prevent reaching limit positions which could damage the mirror or cables.
        
        Args:
            None

        Returns:
            None

        """
        self.timer_count += 1
        self.pos = self.io.getdata()
        
        if self.state == 0:
            self.stop()
            self.real_time_data_window.update_mode("Stop") 
        elif self.state == 1:
            self.home(self.pos)
        elif self.state == 2:
            self.fixed_target() 
        elif self.state == 3:
            self.passive_tracking(self.pos)
        elif self.state == 4:
            self.active_tracking(self.pos)

        elif self.state == 5:
            self.star_tracking(self.pos)
            print("star tracking")
        elif self.state == 6:
            self.manual_control()
        elif self.state == 7:
            self.stop() 
        self.real_time_data_window.update_variables(self.M,self.pos)
        self.t = self.timer_count*self.update_period/3600+self.t0
        self.update_plots()
        self.soft_end_switch = self.io.send_data(self.M)
        self.io.request_motor_data()
        self.logger.log(self.t,self.state,self.pos,self.M,self.E)
        
        
        if self.soft_end_switch == 1:
            self.state = 0
            self.real_time_data_window.update_mode("Reached soft end switch") 
            
            
        
        #Making a limit check to prevent from limit positions 
        if self.pos['A'] < self.low_limit_A and self.M['A'] < 0:
            self.state = 7
        elif self.pos['A'] > self.high_limit_A and self.M['A'] > 0:
            self.state = 7 
            
        if self.pos['E'] < self.low_limit_E and self.M['E'] < 0:
            self.state = 7
        elif self.pos['E'] > self.high_limit_E and self.M['E'] > 0:
            self.state = 7 
            
    def update_plots(self):
        """
        This function updates the plots displayed in different widgets with new data.
        It appends the current time value to the time vector and updates the diagnose plot,
        motion parameters plot, 4Q plot, and spatial view plot with the updated data.

        Args:
            None

        Returns:
            None
        """
        self.t_vector.append(self.t)
        if len(self.t_vector) == 100:
            self.t_vector = self.t_vector[1:]
        self.diagnose.update_diagnose_plot(self.t_vector,self.M)
        self.motion_params.update_motion_plot(self.t_vector, self.pos, self.E)  
        self.q_params.update_fourq_plot(self.t_vector,self.pos)
        self.spatial_view.update_view_plot(self.pos)     
              
    def home(self, pos):
        """
        This function checks if the system is at the home position by comparing the current
        position with the predefined home position with a tolerance of +/- 0.1 degrees for
        both primary and secondary axes. If the system is at the home position, it sets the
        state to 0 (stopped) and updates the real-time data window to indicate that the
        home position is reached and the motors are off. If the system is not at the home
        position, it initiates movement towards the home position using the 'go_to_target'
        function and updates the real-time data window to indicate that the system is
        moving towards the home position.

        Args:
            pos (dict): Current position of the system.

        Returns:
            None
        """
        
        if self.home_pos['E']-0.1 < pos['E'] < self.home_pos['E']+ 0.1 and self.home_pos['A'] - 0.1 < pos['A'] < self.home_pos['A'] + 0.1:
            self.state = 0
            self.real_time_data_window.update_mode("Home position reached, motors off")
        else:
            self.go_to_target(pos, self.home_pos['A'],self.home_pos['E'])
            self.real_time_data_window.update_mode("Moving towards home position")
        
    def fixed_target(self):
        """
        This function calculates the motion required to move the system towards a fixed target position.
        It uses the 'go_to_target' function to initiate movement towards the target position based on the
        current position ('pos') and the target position ('target'). It updates the motion parameters in ('M')
        and updates the real-time data window to indicate that the system is moving towards the target position.

        Args:
            None
            
        Returns:
            dict: The motion parameters required to move towards the fixed target position.
        """
        self.M = self.go_to_target(self.pos, self.target['A'], self.target['E'])
        mode = f"Moving towards {self.target['A']}/{self.target['E']}"
        self.real_time_data_window.update_mode(mode)
        
    def go_to_target(self, POS,TA,TE):
        """
        This function calculates the motion parameters required to move the system from its current position ('POS')
        towards a target position specified by the target angles ('TA' for primary axis and 'TE' for secondary axis).
        It computes the differences between the current position and the target position for both axes ('DA' and 'DE').
        Based on these differences, it determines the motion required for each axis ('M['A']' and 'M['E']').
        If the difference in position is significant (> 2 degrees), the maximum speed ('M_max') is used.
        Otherwise, a slower speed ('M_slow_fixed_target') is calculated proportionally to the difference.
        If the difference is very small (< 0.01 degrees), the motion for that axis is set to 0.
        The calculated motion parameters are returned.

        Args:
            POS (dict): Current position of the system
            TA (float): Target angle for the primary axis.
            TE (float): Target angle for the secondary axis.

        Returns:
            dict: The motion parameters required to move towards the target position.
        """
        DA = POS['A'] - TA
        DE = POS['E'] - TE

        if DE < 0 and DE < -2:
            self.M['E'] = -self.M_max
        elif DE > 0 and DE > 2:
            self.M['E'] = self.M_max
        else:
            self.M['E'] = DE * self.M_slow_fixed_target * self.grad2step

        if DA < 0 and DA < -2:
            self.M['A'] = self.M_max
        elif DA > 0 and DA > 2:
            self.M['A'] = -self.M_max
        elif 0 <= DA < 0.01:
            self.M['A'] = 0
        else:
            self.M['A'] = DA * -self.M_slow_fixed_target * self.grad2step
        return self.M
    
    def passive_tracking(self, pos):
        """  
        This function implements passive tracking, which adjusts the system's motion based on the difference
        between the current position and a reference position (position that directs the beam of the sun into the laboratory). 
        It calculates the differences ('DA' and 'DE') between the current position ('pos') and the reference position adjusted with an offset.
        It then adjusts the system's motion using proportional-integral (PI) controllers for both azimuth and
        elevation axes, taking into account the passive speed settings and the current state of the system.
        The system adjusts its motion differently based on whether the differences are significant or within
        fine-tuning thresholds. It updates the system's state and motion parameters accordingly and updates
        the real-time data window to indicate the tracking mode.

        Args:
            pos (dict): Current position of the system.

        Returns:
            None
        """
        self.pos_and_speed()
        DA = pos['A'] - (self.pos1['A']-self.target_offset['A'])
        DE = pos['E'] - (90+(self.pos1['E']-self.target_offset['E']))/2

        self.passive_speed['A'] = 2*self.passive_speed['A']
        self.passivPI_E.feed(DA, self.passive_speed['A'])
        self.passivPI_A.feed(DE, self.passive_speed['E'])
        
        #Elevation
        if DE < 0 and DE < -1:
            self.M['E'] = -self.M_max
            self.state_secondary = 1
        elif DE > 0 and DE > 1:
            self.M['E'] = self.M_max
            self.state_secondary = 1
        elif DE > 0.01 or DE < -0.01:
            self.state_secondary = 0
            self.M['E'] = -(self.passive_speed['E']*self.grad2step-DE*self.grad2step*self.M_slow_rel)
        else:
            self.state_secondary = 0
            self.M['E'] = -self.passivPI_E.get_speed()
        
        #Azimuth   
        if DA < 0 and DA < -1:
            self.state_primary = 1
            self.M['A'] = self.M_max
        elif DA > 0 and DA > 1:
            self.state_primary = 1
            self.M['A'] = -self.M_max 
        elif DA > 0.01 or DA <-0.01:
            self.state_primary = 0
            self.M['A'] = self.passive_speed['A']*self.grad2step-DA*self.grad2step*self.M_slow_rel_Az
        else:
            self.state_primary = 0
            self.M['A'] = self.passivPI_A.get_speed()
        
        if self.state_secondary == 1 and self.state_primary == 1:
            self.real_time_data_window.update_mode("Passive tracking: coarse")
            
        elif self.state_secondary == 0 and self.state_primary == 0:
            self.real_time_data_window.update_mode("Passive tracking: fine")
            
        
        self.E['DA'] = DA
        self.E['DE'] = DE
    
    def active_tracking(self, pos):
        """
        This function implements active tracking, which adjusts the system's motion based on feedback from
        the 4Q sensor. It calculates the differences ('DA' and 'DE') between the current position ('pos')
        and a reference position adjusted with an offset. It computes the orientation angles ('DAO' and 'DEO')
        from the 4Q sensor feedback and adjusts the system's motion using proportional-integral (PI)
        controllers for both azimuth and elevation axes, taking into account the active sensor feedback
        and the system's state. The function evaluates whether the system is within the active tracking frame
        and updates the tracking mode accordingly. If the system is out of frame, it switches to passive tracking.
        It then adjusts the system's motion based on the calculated speeds and updates the motion parameters ('M').
        Finally, it updates the error values for logging purposes.

        Args:
            pos (dict): Current position of the system.
            
        Returns:
            None
        """
        
        #print(pos['FQ'][0])
        #print(pos['FQ'][1])
        
        DAO = np.arctan(pos['FQ'][0])/ np.sin(np.deg2rad(90-self.pos['E']))/2
        DEO = np.arctan(pos['FQ'][1])/2
        
        if DAO == np.nan:
            DAO = 0
        if DEO == np.nan:
            DEO = 0
        
        self.pos_and_speed()
        DA = pos['A'] - (self.pos1['A']-self.target_offset['A'])
        DE = pos['E'] - (90+(self.pos1['E']-self.target_offset['E']))/2
        
        self.passive_speed['A'] = 2*self.passive_speed['A']
        
        self.activePI_A.feed(DAO, self.passive_speed['A'])
        self.activePI_E.feed(DEO, self.passive_speed['E'])
        
        if abs(DE) < self.active_frame and abs(DA) < self.active_frame:
            self.active_tracking_condition = True
        else:
            self.active_tracking_condition = False
        

        if self.active_tracking_condition:
            if np.all(self.activePI_E.array == 0):
                self.real_time_data_window.update_mode(f"No active signal, tracking still in frame {self.active_frame} deg")
            else:
                self.real_time_data_window.update_mode(f"Active tracking in frame, frame size {self.active_frame} deg")    
        else:
            self.real_time_data_window.update_mode("Out of frame: switch to passive tracking")
            #self.state = 3 # Fraglich ob überhaupt benötigt wird weil automatischer wechsel zwischen active und passive tracking stattfindet in der active tracking funktion
        
        if DEO < 0 and DEO < -2:
            self.M['E'] = self.M_max
        elif DEO > 0 and DEO > 2:
            self.M['E'] = -self.M_max
        elif self.active_tracking_condition == False:
            self.M['E'] = -(self.passive_speed['E']*self.grad2step-DE*self.grad2step*self.M_slow_rel)
        else:
            self.M['E'] = -self.activePI_E.get_speed()
            
        if DAO < 0 and DAO < -2:
            self.M['A'] = -self.M_max
        elif DAO > 0 and DAO > 2:
            self.M['A'] = self.M_max
        elif self.active_tracking_condition == False:
            self.M['A'] = self.passive_speed['A']*self.grad2step-DA*self.grad2step*self.M_slow_rel_Az
        elif abs(DA) > 0.2:
            self.M['A'] = self.passive_speed['A']*self.grad2step+DAO*self.grad2step*self.M_slow_rel_Az
        else:
            self.M['A'] = self.activePI_A.get_speed()
            
        self.E['DA'] = DA
        self.E['DE'] = DE
        
    def star_tracking(self, pos):
        """
        This function implements star tracking, which adjusts the system's motion to track a celestial object
        based on its right ascension (RA) and declination (Dec). It calculates the azimuth and zenith
        angles of the star using the target RA and Dec and transforms them into the system's coordinate system.
        The function then calculates the differences ('DA' and 'DE') between the current position ('pos') and
        the calculated position of the star. Based on these differences, it adjusts the system's motion to
        track the star using proportional control. The motion parameters ('M') are updated accordingly.

        Args:
            pos (dict): Current position of the system.

        Returns:
            None
        """
        azimuth, zenith = self.get_angles_star(self.target_ra, self.target_dec)
        star = {'A': azimuth, 'E':90-zenith}
        Heli_star = PRM.polar_rotation(star['A'],90-star['E'], self.trafoA, self.trafoE)
        DA = pos['A'] - Heli_star['A']
        DE = pos['E'] - (90+Heli_star['E'])/2
        
        
        #Fraglich ob star tracking so funktioniert oder besser gleiches prinzip wie im passive tracking
        if DE < 0 and DE < -1:
            self.M['E'] = -self.M_max   
        elif DE > 0 and DE > 1:
            self.M['E'] = self.M_max
        else:
            self.M['E'] = DE * -self.M_slow_rel
            
        if DA < 0 and DA < -1:
            self.M['E'] = -self.M_max   
        elif DA > 0 and DA > 1:
            self.M['E'] = self.M_max
        else:
            self.M['E'] = DA * -self.M_slow_rel
            
            
    def manual_control(self):
        """
        This function enables manual control of the system based on user input. It requests manual input
        from the input/output interface and adjusts the system's motion accordingly.
        It checks the status of each motor and adjusts the system's motion speed based on the received commands.
        The motion parameters ('M') are updated accordingly for both azimuth ('A') and elevation ('E') axes.

        Args:
            None
            
        Returns:
            None
        """
        self.io.request_manual_input()
        #print(self.io.M_rep_M)
        try:
            statusE = bin(int(self.io.M_rep_M['E']))
            if statusE[-1] == "1":
                #print("ror")
                self.M['E'] = -self.speed
            elif statusE[-2] == "1":
                #print("rol")
                self.M['E'] = self.speed
            else:
                #print("stop")
                self.M['E'] = 0
        except:
            print("No valid reply from motor E")
            self.M['E'] = 0
        
        try:
            statusA = bin(int(self.io.M_rep_M['A']))
            if statusA[-1] == "1":
                #print("ror")
                self.M['A'] = -self.speed
            elif statusA[-2] == "1":
                #print("rol")
                self.M['A'] = self.speed
            else:
                #print("stop")
                self.M['A'] = 0
        except:
            print("No valid reply from motor A")
            self.M['A'] = 0
        
    def pos_and_speed(self):
        """
        This function calculates the current position and speed of the system based on the current time.
        It obtains the azimuth and zenith angles at two consecutive time points and converts them into
        the system's coordinate system using polar rotation. It then calculates the speed of the system
        by taking the difference in position between the two time points which is further used in passive tracking. 
        
        Args:
            None
            
        Returns:
            None
        """

        now = datetime.now()
        azimuth, zenith = self.get_angles(now)
        self.pos1 = PRM.polar_rotation(self, azimuth,90-zenith, self.trafoA, self.trafoE)
        
        now1 = datetime.now() + timedelta(seconds=1)
        azimuth, zenith = self.get_angles(now1)
        self.pos2 = PRM.polar_rotation(self, azimuth,90-zenith, self.trafoA, self.trafoE)
    
        self.passive_speed['A'] =  self.pos2['A'] - self.pos1['A']   # speed calculation in grad/second (important for PI Controller afterwards)
        self.passive_speed['E'] =  self.pos2['E'] - self.pos1['E']   # speed calculation in grad/second (important for PI Controller afterwards)
        
    def get_angles(self, time):
        """
        This function calculates the solar azimuth and zenith angles based on the given time.
        It uses the 'get_solarposition' function to obtain the solar position for the specified time
        and location. The azimuth and zenith angles are extracted from the solar position data.

        Args:
            time (datetime): The time at which the solar angles are to be calculated.

        Returns:
            tuple: A tuple containing the solar azimuth and zenith angles.
        """
        timestamp = pd.Timestamp(time, tz='Europe/Zurich')
        solar_position = get_solarposition(timestamp, self.latitude, self.longitude, self.altitude)
        azimuth = solar_position['azimuth'].values
        zenith = solar_position['apparent_zenith'].values
        return azimuth,zenith
    
    def get_angles_star(self, target_ra, target_dec):
        """
        This function calculates the azimuth and zenith angles for a celestial object based on
        its right ascension (RA) and declination (Dec). It creates a SkyCoord object for the specified
        RA and Dec, transforms it to an AltAz frame, and obtains the azimuth and altitude angles.
        The calculated azimuth and zenith angles are returned.

        Args:
            target_ra (float): Right ascension of the target celestial object (in degrees).
            target_dec (float): Declination of the target celestial object (in degrees).

        Returns:
            tuple: A tuple containing the azimuth and zenith angles of the celestial object.
        """
        position = SkyCoord(ra=target_ra*u.degree, dec=target_dec*u.degree, frame='icrs')
        location = EarthLocation(lat=self.latitude*u.deg, lon=self.longitude*u.deg, height=self.altitude*u.m)
        time = Time(datetime.now())

        altaz = position.transform_to(AltAz(obstime=time,location=location))
        alt = altaz.alt.deg
        az = altaz.az.deg
        return az, alt
    
    def stop(self):
        """
        This function stops the motion of the system.

        Args:
            None
            
        Returns:
            None
        """
        self.M['E'] = 0
        self.M['A'] = 0
        
        
       
    #Things to do when the application gets closed   
    def closeEvent(self,event):
        """
        This function handls the close event of the main window.

        Args:
            event: The close event triggered by closing the main window.

        Returns:
            None
        """
        self.real_time_data_window.close()
        self.motion_params.close()
        self.q_params.close()
        self.diagnose.close()
        self.spatial_view.close()
        self.timer.stop()
        self.timer_fourq.stop()
        self.fourq.delete()
        self.io.delete()
        event.accept()
        
        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = HeliostatUI()
    window.show()
    sys.exit(app.exec_())