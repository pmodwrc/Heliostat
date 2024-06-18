from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import os
import numpy as np
from sympy import Point3D, Line3D



class PRM(QObject):
    """
    This class calculates the vectors and lines for the 2D/3D view of the Heliostat.
    
    Attributes:
        nn (int): Number of drawing points for a circle.
        m1r (int): Radius of the primary mirror
        m2r (int): Radius of the secondary mirror.
        AU (float): Sun-earth distance
        trafo_A_haus (float): Angle between Heliostat and house normal vector
        trafoA (float): Coordinate transformation parameter AzEl2Heli.
        trafoE (float): Coordinate transformation parameter AzEl2Heli.
        n_p1 (numpy.ndarray): Normal vector of the primary mirror.
        PV (numpy.ndarray): View point in the laboratory.
        N12 (numpy.ndarray): 
        d2l (int): 
        d12 (int): 
        P2c (numpy.ndarray): Center of secondary mirror.
        P1c (numpy.ndarray): Center of primary mirror.
        n_p2 (numpy.ndarray): Normalvector of mirror 2
        n_pl (numpy.ndarray): Normalvector of the projection plane.
    """
    def __init__(self):
        """
        Initialize attributes.
        
        Retunrs:
            None
        """
        self.nn = 101
        self.m1r = 315
        self.m2r = 300
        self.AU = 1.4960e14
        
        self.trafo_A_haus = -61.1892
        self.trafoA = 4.91
        self.trafoE = 26.9
        
        self.n_p1 = np.array([0.672503844142359, -0.233550637023614, 0.702276782728587])
        self.PV = np.array([8700,1500,1500])
        
        self.N12 = self.drehen(self.trafoE,[-1,0,0],[0,0,1])
        self.N12 = self.drehen(self.trafo_A_haus,[0,0,1],self.N12)
        
        self.d2l = 2500
        self.d12 = 1010
        self.P2c = self.PV+[0, -self.d2l,0]
        self.P1c = self.P2c- self.N12*self.d12
        
        self.n_p2 = ((self.P1c-self.P2c)/np.linalg.norm(self.P1c-self.P2c)) + ((self.PV-self.P2c)/np.linalg.norm(self.PV-self.P2c))
        self.n_pl = np.array([0,1,0])
        
    def update_position(self, pos):
        """
        This method updates the position of the PRM object based on the azimuth ('A') and 
        elevation ('E') angles provided in the 'pos' dictionary. It calculates the normal 
        vector 'n_p1' corresponding to mirror 1, as well as the direction of 
        the sun 'V_sun' and the sun position 'Pc_sun'.

        Args:
            pos (dict): A dictionary containing 'A' and 'E' keys representing azimuth 
                        and elevation angles, respectively.

        Returns:
            None
        """
        self.heli = self.inverse_polar_rotation(pos['A'], pos['E'], -self.trafo_A_haus, self.trafoE)
        A_rad= np.deg2rad(self.heli['A'])
        E_rad= np.deg2rad(self.heli['E']) 
        self.n_p1[0], self.n_p1[1], self.n_p1[2] = self.sph2cart(A_rad, E_rad, 1)
        
        self.V_sun = self.P1c-2*self.intersection_point(self.P2c,self.n_p1,self.P1c,self.n_p1) + self.P2c
        self.V_sun = self.V_sun/np.linalg.norm(self.V_sun)
        self.Pc_sun = self.P1c + self.V_sun*self.AU

    def intersection_point(self,P1, v, P2, n):
        """
        Calculate the intersection point of a line and a plane.

        Args:
            P1 (numpy.ndarray): Coordinates of a point on the line.
            v (numpy.ndarray): Direction vector of the line.
            P2 (numpy.ndarray): Coordinates of a point on the plane.
            n (numpy.ndarray): Normal vector of the plane.

        Returns:
            numpy.ndarray: Coordinates of the intersection point.
        """
        
        n = -1*n
        D = np.dot(n, P2)
        t = (D - np.dot(n, P1)) / np.dot(n, v)
        s = P1 + t * np.array(v)
        return s  
    
    @staticmethod
    def polar_rotation(self, A,E,alpha,beta):
        """
        This method applies polar rotation to given azimuth and elevation angles (A, E) 
        by rotating them using the angles alpha and beta. The rotation is performed 
        around the z-axis followed by the y-axis. 

        Args:
            A (float): Azimuth angle in degrees.
            E (float): Elevation angle in degrees.
            alpha (float): Rotation angle around the z-axis in degrees.
            beta (float): Rotation angle around the y-axis in degrees.

        Returns:
            dict: A dictionary containing the updated azimuth and elevation angles 
                after the polar rotation.
        """

        A_rad = np.deg2rad(A)
        E_rad = np.deg2rad(E) 
        x, y, z = PRM.sph2cart(PRM,A_rad, E_rad,1)
        v_out = PRM.drehen(PRM,alpha,np.array([0,0,1]),np.array([x[0],y[0],z[0]]))
        v_out = PRM.drehen(PRM,beta, np.array([0,-1,0]),v_out)
        x, y, z = v_out

        A_new_rad, E_new_rad = PRM.cart2sph(PRM,x, y, z)
        A_new = np.degrees(A_new_rad)
        E_new = np.degrees(E_new_rad)

        if A_new < 0:
            A_new += 360
        Heli = {'A': A_new, 'E': E_new}
        return Heli  
       
    def inverse_polar_rotation(self, A,E,alpha,beta):
        """
        This method applies inverse polar rotation to given azimuth and elevation angles (A, E) 
        by first rotating them around the z-axis by alpha degrees followed by a rotation around 
        the x-axis by beta degrees.

        Args:
            A (float): Azimuth angle in degrees.
            E (float): Elevation angle in degrees.
            alpha (float): Rotation angle around the z-axis in degrees.
            beta (float): Rotation angle around the x-axis in degrees.

        Returns:
            dict: A dictionary containing the updated azimuth and elevation angles 
                after the inverse polar rotation.
        """
        A = 270 - A
        E= 180 - E
        A_rad = np.deg2rad(A)
        E_rad = np.deg2rad(E) 
        x, y, z = self.sph2cart(A_rad, E_rad,1)
        v_out = self.drehen(-beta,np.array([1,0,0]),np.array([x,y,z]))
        v_out = self.drehen(-alpha, [0,0,1],v_out)
        x, y, z = v_out

        A_new_rad, E_new_rad = self.cart2sph(x,y,z)
        A_new = np.degrees(A_new_rad)
        E_new = np.degrees(E_new_rad)
    
        if A_new < 0:
            A_new += 360
        Heli = {'A': A_new, 'E': E_new}
        return Heli  
    
    def cart2sph(self,x,y,z):
        """
        Convert cartesian coordinates to spherical coordinates.

        This method calculates the azimuth, elevation, and radial distance (r) 
        from the origin to a point defined by its cartesian coordinates (x, y, z).

        Args:
            x (float): X-coordinate of the point.
            y (float): Y-coordinate of the point.
            z (float): Z-coordinate of the point.

        Returns:
            tuple: A tuple containing the azimuth angle (in radians), 
                elevation angle (in radians), and radial distance (r).
        """
        azimuth = np.arctan2(y,x)
        elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    
        return azimuth, elevation
        
    def sph2cart(self, azimuth,elevation,r):
        """
        Convert spherical coordinates to cartesian coordinates.
        
        This method calculates the cartesian coordinates (x, y, z) 
        from given spherical coordinates: azimuth angle, elevation angle, and radial distance (r).

        Args:
            azimuth (float): Azimuth angle in radians.
            elevation (float): Elevation angle in radians.
            r (float): Radial distance from the origin to the point.

        Returns:
            tuple: A tuple containing the x, y, and z coordinates of the point in cartesian space.
        """
        x = r * np.cos(elevation) * np.cos(azimuth)
        y = r * np.cos(elevation) * np.sin(azimuth)
        z = r * np.sin(elevation)
        return x, y, z
    
    def project_m1(self, pos):
        """
        Project mirror M1 onto the same plane as mirror M2.

        This method updates the position of mirror M1 based on the given position 'pos', 
        then projects M1 onto the same plane as M2. 

        Args:
            pos (dict): A dictionary containing 'A' and 'E' keys representing azimuth and elevation angles.

        Returns:
            numpy.ndarray: An array containing the coordinates of the projected points of mirror M1.
        """
        self.update_position(pos)
        msp1 = self.punkt_spiegeln_einfach(self.P1c)
        nsp1 = self.vektor_spiegeln_einfach(self.n_p1,self.P1c)
        Abb = self.Abbildung(self.circle_points(msp1, nsp1, self.m1r, self. nn))
        return Abb
    def project_m2(self):
        """
        Project mirror M2 onto its own plane.

        This method projects mirror M2 onto the wall of the labratory.

        Returns:
            numpy.ndarray: An array containing the coordinates of the projected points of mirror M2.
        """
        Abb = self.Abbildung(self.circle_points(self.P2c,self.n_p2, self.m2r, self.nn))
        return Abb
    
    
    
    def Abbildung(self, points):
        """
        This function calculates the orthogonal projection of a given object in 3D space
        """
        num_rows, num_collums = np.shape(points)
        AbbZ = np.zeros_like(points)
        for i in range(0,num_rows):
            AbbZ[i] = self.intersection_point(self.PV, self.PV-points[i],self.PV-self.n_pl,self.n_pl)-self.PV       
        return AbbZ
     
    def punkt_spiegeln_einfach(self,P_in):   
        """

        This method mirrors the given point 'P_in' on the plane defined by mirror M2.

        Args:
            P_in (numpy.ndarray): Coordinates of the point to be reflected.

        Returns:
            numpy.ndarray: Coordinates of the reflected point.
        """
        P_out = self.intersection_point(P_in, self.n_p2, self.P2c, self.n_p2)
        P_out = P_out + (P_out - P_in)
        return P_out
    
    def vektor_spiegeln_einfach(self, V_in,  P_in ):
        """
        This method mirrrors the given vector 'V_in' on the plane defined by mirror M2.

        Args:
            V_in (numpy.ndarray): Vector to be reflected.
            P_in (numpy.ndarray): Coordinates of the point on which the reflection is performed.

        Returns:
            numpy.ndarray: Coordinates of the reflected vector.
        """
        P_m = self.intersection_point(P_in, self.n_p2, self.P2c, self.n_p2)
        P_m = P_m + (P_m-P_in)
        V_out = self.intersection_point(P_in,V_in, self.P2c, self.n_p2)
        V_out = -(V_out-P_m)
        return V_out     
     
       
    def circle_points(self, P,N,r,n):
        """
        Calculate points on a circle with given center, radius, and normal vector.

        This method calculates 'n' points on a circle with center 'P', radius 'r', and normal vector 'N'.

        Args:
            P (numpy.ndarray): Coordinates of the center of the circle.
            N (numpy.ndarray): Normal vector of the plane containing the circle.
            r (float): Radius of the circle.
            n (int): Number of points to generate on the circle.

        Returns:
            numpy.ndarray: An array containing the coordinates of the points on the circle.
        """
        N = N/np.linalg.norm(N)
        NR = np.cross(N,[1312, 213512,3415]) 
        NR = NR/np.linalg.norm(NR)*r
        CP = np.empty((n,3))
        for i in range(0,n):
            CP[i]= np.array(P) + self.drehen(360/(n)*i,N,NR)
        return CP
        
    def drehen(self,alpha,n,v_in):
        """
        Rotate a vector around a given axis with a specified angle.

        This method rotates the given vector 'v_in' around the axis 'n' by the angle 'alpha'.

        Args:
            alpha (float): Angle of rotation in degrees.
            n (numpy.ndarray): Axis of rotation.
            v_in (numpy.ndarray): Vector to be rotated.

        Returns:
            numpy.ndarray: The rotated vector.
        """
        n = np.array(n)/np.linalg.norm(np.array(n))
        v1 = np.dot(n, np.dot(np.transpose(n), v_in)) 
        v2 = np.cross(np.cross(n, v_in), n) *np.round(np.cos(np.deg2rad(alpha)),5)
        v3 = np.cross(n, v_in) * np.round(np.sin(np.deg2rad(alpha)),5)
        v_out = v1+v2+v3
        return v_out
        
        
        
    
    
    
    