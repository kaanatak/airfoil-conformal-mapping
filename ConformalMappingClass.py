import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from ConformalMappingFunctions import *

class ConformalMapping:
    def __init__(self, airfoil_data=None, xc=None, yc=None, radius=None,U_inf=1, aoa=20, c=1, N_r=100, N_theta=145):
        """
        Initialize the ConformalMapping object. You can either pass airfoil coordinates or circle parameters.
        
        :param airfoil_coords: Airfoil coordinates as complex numbers (optional)
        :param xc: Circle center x-coordinate (optional)
        :param yc: Circle center y-coordinate (optional)
        :param c: Joukowski transformation constant (optional)
        :param radius: Radius of the circle (optional)
        :param N_r: Number of radial points in the grid (default 100)
        :param N_theta: Number of angular points in the grid (default 100)
        """

        # Grid parameters
        self.N_r = N_r
        self.N_theta = N_theta

        # Joukowski transformation parameters
        self.c = c

        # Flow parameters
        self.U_inf = U_inf
        self.AoA = aoa*np.pi/180
        
        if airfoil_data is not None:
            
            airfoil_coords = airfoil_data[:, 0] + 1j * airfoil_data[:, 1]

            initial_guess = [0.0, 0.0, 1.0]  
            # Joukowski parameters
            c = 1

            result = minimize(objective, initial_guess, args=(airfoil_coords,c))
                
            x_c_opt, y_c_opt, R_opt = result.x
            self.xc = x_c_opt
            self.yc = y_c_opt
            self.radius = R_opt
            self.airfoil_coords = airfoil_coords
            print(f"Optimized Circle Parameters: x_c = {x_c_opt}, y_c = {y_c_opt}, R = {R_opt}")
        else:
            # If airfoil_coords is not provided, use the provided circle parameters
            self.xc = xc
            self.yc = yc
            self.radius = radius
            self.airfoil_coords = None
        print(f"Circle Parameters: x_c = {self.xc}, y_c = {self.yc}, R = {self.radius}")
        print(f"Joukowski Transformation Constant: c = {self.c}")
        print(f"Flow Parameters: U_inf = {self.U_inf}, AoA = {self.AoA}")
        print(f"Grid Parameters: N_r = {self.N_r}, N_theta = {self.N_theta}")

        # Initialize the variables to store the results
        self.results = {}

    
    def solve(self):

        # Circle parameters
        xc = self.xc
        yc = self.yc
        a =  self.radius

        # Grid parameters
        N_r = self.N_r
        N_theta = self.N_theta
        
        # Joukowski Transformation parameters
        c = self.c

        # Flow parameters
        AoA = self.AoA
        U_inf = self.U_inf

        circle, airfoil = geometry(xc,yc,a,N_theta,c)

        Z = generate_grid(N_r,N_theta,xc,yc,a)

        Z_transformed = joukowski_transformation(Z,c)

        Kappa = 2*np.pi*U_inf*a**2
        Gamma = 4*np.pi*U_inf*a*np.sin(AoA)

        # Doublet Coordinates
        xd = xc
        yd = yc

        # Vortex Coordinates
        xv = xc
        yv = yc

        Phi = uniform_potential(Z,U_inf) + doublet_potential(Z, Kappa, xd, yd) + vortex_potential(Z, Gamma, xv, yv)

        u_uni, v_uni = uniform_flow(Z, U_inf, AoA)
        u_d, v_d = doublet_flow(Z, Kappa, xc, yc)
        u_v, v_v = vortex_flow(Z, Gamma, xc, yc)

        u_z = u_uni + u_d + u_v
        v_z = v_uni + v_d + v_v

        u_transformed, v_transformed = transform_velocity(u_z,v_z,Z,c)

        Cp = compute_pressure(u_z,v_z,U_inf)
        Cp_transformed = compute_pressure(u_transformed,v_transformed,U_inf)

        L_circle = calculate_lift(Z, Cp,Z,yc)
        L_airfoil = calculate_lift(Z_transformed, Cp_transformed,Z,yc)

        self.results['Z'] = Z
        self.results['circle'] = circle
        self.results['Z_transformed'] = Z_transformed
        self.results['airfoil'] = airfoil
        self.results['Phi'] = Phi
        self.results['Cp'] = Cp
        self.results['Cp_transformed'] = Cp_transformed
        self.results['L_circle'] = L_circle
        self.results['L_airfoil'] = L_airfoil
        self.results['AoA'] = AoA
        self.results['u_z'] = u_z
        self.results['v_z'] = v_z
        self.results['u_transformed'] = u_transformed
        self.results['v_transformed'] = v_transformed


    def get_results(self):
        return self.results
    
    def plot(self, features):

        if isinstance(features, str):
            features = [features]

        for feature in features:
            if feature == 'potential':
                plot_potential(self.xc, self.yc, self.radius, self.results['Z'], self.results['Phi'], self.results['circle'],self.results['Z_transformed'], self.results['airfoil'],self.c,self.AoA)
            elif feature == 'velocity':
                plot_velocity(self.results['Z'], self.results['u_z'], self.results['v_z'], self.results['circle'], self.results['Z_transformed'], self.results['u_transformed'], self.results['v_transformed'], self.results['airfoil'],stride=5)
            elif feature == 'pressure':
                plot_pressure(self.results['Z'], self.results['Cp'], self.results['circle'], self.results['Z_transformed'], self.results['Cp_transformed'], self.results['airfoil'])
            elif feature == 'pressure_distribution':
                plot_pressure_difference(self.results['Z'], self.results['Cp'], self.results['Z_transformed'], self.results['Cp_transformed'],self.yc)
            elif feature == 'grid':
                plot_grid(self.results['Z'], self.results['circle'], self.results['Z_transformed'], self.results['airfoil'])
            elif feature == 'fitting':
                if self.airfoil_coords is not None:
                    plot_fitting(self.airfoil_coords, self.results['airfoil'])
                else:
                    print("Airfoil coordinates are not provided. Cannot plot fitting.")
            elif feature == 'general':
                plot_all(self.results['Z'], self.results['circle'], self.results['Z_transformed'], self.results['airfoil'], self.results['Phi'], self.results['Cp'], self.results['Cp_transformed'], self.yc, stride=5, AoA=self.AoA)        
            elif feature == 'mapping':
                plot_mapping(self.results['circle'], self.results['airfoil'])
            else:
                print("Invalid feature. Available options: 'grid', 'fitting', 'potential', 'velocity', 'pressure', 'pressure_distribution', 'general'")
                print("You can also pass a list of features to plot multiple plots.")
        plt.show()

    def interactive(self):
        interactive_plot()
