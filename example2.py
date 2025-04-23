import ConformalMappingClass as cm
import numpy as np

# Initialize ConformalMapping class with some values (values not important) for interactive plot
x = -0.15
y = 0.23
r = 0.23*np.sqrt(13*2)

# Create an instance of the ConformalMapping class
cmapping = cm.ConformalMapping(xc=x, yc=y, radius=r,N_r=200, N_theta=250)

# Call the interactive method to plot the interactive plot
cmapping.interactive()