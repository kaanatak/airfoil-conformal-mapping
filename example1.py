import ConformalMappingClass as cm
import numpy as np

# Load airfoil data from file (x,y) coordinates
airfoil_coord_data = np.loadtxt("airfoil_data2.csv", delimiter=",")

# Initialize ConformalMapping class with airfoil data --> this will optimize the airfoil data to get circle parameters
cmapping = cm.ConformalMapping(airfoil_data=airfoil_coord_data)

# Solve the problem
cmapping.solve()

# Plot the fitting after optimization
cmapping.plot('fitting')

# Plot the solutions
cmapping.plot('mapping')

# Plot the streamlines
cmapping.plot('streamlines')

# Plot the velocity field
cmapping.plot('velocity')

# Plot the pressure coefficient
cmapping.plot('pressure')

# Plot pressure distribution on the airfoil
cmapping.plot('pressure_distribution')

