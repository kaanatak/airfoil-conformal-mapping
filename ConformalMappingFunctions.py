import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox, RadioButtons

def joukowski_transformation(z,c):
    return z + c**2 / z

def d_joukowski_dz(z, c):
    return 1 - c**2 / z**2

def uniform_potential(Z, U_inf):
    return U_inf*Z.imag

def doublet_potential(Z, Kappa, xc, yc):
    return - Kappa / (2 * np.pi) * (Z.imag-yc) / ((Z.real-xc)**2 + (Z.imag - yc)**2)

def vortex_potential(Z, Gamma, xc, yc):
    return Gamma  / (4 * np.pi) * np.log((Z.real - xc)**2 + (Z.imag - yc)**2)

def uniform_flow(Z, U_inf, AoA):
    u_z = U_inf*np.cos(AoA)
    v_z = U_inf*np.sin(AoA)
    return u_z, v_z

def doublet_flow(Z, Kappa, xd, yd):
    u_z = (- Kappa / (2 * np.pi) * ((Z.real - xd)**2 - (Z.imag - yd)**2) / ((Z.real - xd)**2 + (Z.imag - yd)**2)**2)
    v_z = (- Kappa/ (2 * np.pi) * 2 * (Z.real - xd) * (Z.imag - yd) / ((Z.real - xd)**2 + (Z.imag - yd)**2)**2)
    return u_z, v_z

def vortex_flow(Z, Gamma, xv, yv):
    u_z = Gamma / (2 * np.pi) * (Z.imag - yv) / ((Z.real - xv)**2 + (Z.imag - yv)**2)
    v_z = - Gamma / (2 * np.pi) * (Z.real - xv) / ((Z.real - xv)**2 + (Z.imag - yv)**2)
    return u_z, v_z

def generate_grid(N_r,N_theta,xc,yc,a):
    r = np.linspace(a,6,N_r)
    theta = np.linspace(0,2*np.pi,N_theta)
    R,Theta = np.meshgrid(r,theta)
    Z = (xc+1j*yc) + R*np.exp(1j*Theta)
    return Z

def rotate_grid(Z,xc,yc,AoA):
    Z_rotated = -(Z.real-xc)*np.cos(AoA) - (Z.imag-yc)*np.sin(AoA) + 1j*(-(Z.real-xc)*np.sin(AoA) + (Z.imag-yc)*np.cos(AoA)) + xc + 1j*yc
    return Z_rotated

def geometry(xc,yc,a,N_theta,c):
    theta = np.linspace(0,2*np.pi,N_theta)
    circle = (xc+1j*yc) + a*np.exp(1j*theta)
    airfoil = joukowski_transformation(circle,c)
    return circle, airfoil

def transform_velocity(u_z,v_z,Z,c):
    W_z = u_z - 1j*v_z
    W_transformed = W_z/d_joukowski_dz(Z,c)

    u_transformed = W_transformed.real
    v_transformed = -1*W_transformed.imag
    return u_transformed, v_transformed

def compute_pressure(u_z,v_z,U_inf):
    Cp = 1 - (u_z**2 + v_z**2) / U_inf**2
    return Cp

def airfoil_split(Z_transformed,Cp_transformed,Z,yc):
    # Upper Surface
    upper_indices = np.where(Z[:,0].imag >= yc)
    lower_indices = np.where(Z[:,0].imag < yc)
    X_upper = Z_transformed.real[upper_indices,0]
    X_lower = Z_transformed.real[lower_indices,0]
    Cp_upper = Cp_transformed[upper_indices,0]
    Cp_lower = Cp_transformed[lower_indices,0]
    return X_upper, X_lower, Cp_upper, Cp_lower

def calculate_lift(Z_transformed, Cp_transformed,Z,yc):
    X_upper, X_lower, Cp_upper, Cp_lower = airfoil_split(Z_transformed,Cp_transformed,Z,yc)
    
    # mask inf nan values
    X_lower  = X_lower[:,1:-2]
    X_upper = X_upper[:,1:-2]
    Cp_upper = Cp_upper[:,1:-2]
    Cp_lower = Cp_lower[:,1:-2]

    area_upper = np.trapz(Cp_upper,X_upper)
    area_lower = np.trapz(Cp_lower,X_lower)

    L = area_upper - area_lower
    L = L[0]
    return L

def objective(params, airfoil_coords,c):
    x_c, y_c, R = params  # Circle center and radius
    z_c = x_c + 1j * y_c
    
    # Generate points on the circle
    theta = np.linspace(0, 2 * np.pi, len(airfoil_coords))
    circle_points = z_c + R * np.exp(1j * theta)
    
    # Apply Joukowski transform
    airfoil_fit = joukowski_transformation(circle_points,c)
    
    # Compute error between transformed circle and airfoil data
    error = np.sum(np.abs(airfoil_fit - airfoil_coords)**2)
    return error

def plot_potential(xc,yc,a,Z, Phi, circle, Z_transformed, airfoil, c=1 ,AoA=0):

    if AoA != 0:
        ZZ = rotate_grid(Z,xc,yc,AoA)
        ZZ_transformed = joukowski_transformation(ZZ,c)  
    else:
        ZZ = Z
        ZZ_transformed = Z_transformed
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    axs[0].contour(ZZ.real, ZZ.imag, Phi, levels=50)
    axs[0].plot(np.real(circle), np.imag(circle), color='red')
    axs[0].set_title('Original Potential')
    

    axs[1].contour(ZZ_transformed.real, ZZ_transformed.imag, Phi, levels=50)
    axs[1].plot(np.real(airfoil), np.imag(airfoil), color='red')
    axs[1].set_title('Transformed Potential')
    axs[1].axis('equal')
    plt.tight_layout()



def plot_mapping(circle, airfoil):

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    axs[0].plot(np.real(circle), np.imag(circle), color='red')
    #axs[0].set_title('Original Domain')
    axs[0].tick_params(axis='both', labelsize=16)  # Increase tick font size

    axs[1].plot(np.real(airfoil), np.imag(airfoil), color='red')
    #axs[1].set_title('Mapped Domain')
    axs[1].axis('equal')
    axs[1].tick_params(axis='both', labelsize=16)  # Increase tick font size

    plt.tight_layout()

def plot_velocity(Z, u_z, v_z, circle, Z_transformed, u_transformed, v_transformed, airfoil, stride=5):

        s1 = 'default'
        s2 = 20
        yl = -3
        yu = 3
        xl = -3
        xu = 3

        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        axs[0].quiver(Z.real[::stride, ::stride], Z.imag[::stride, ::stride], u_z[::stride, ::stride], v_z[::stride, ::stride])
        axs[0].plot(np.real(circle), np.imag(circle), color='red')
        axs[0].set_title('Original Velocity')
        axs[0].axis('equal')
        axs[0].set_xlim(xl,xu)
        axs[0].set_ylim(yl,yu)
    
        axs[1].quiver(Z_transformed.real[::stride, ::stride], Z_transformed.imag[::stride, ::stride], u_transformed[::stride, ::stride], v_transformed[::stride, ::stride], scale=s2)
        axs[1].plot(np.real(airfoil), np.imag(airfoil), color='red')
        axs[1].set_title('Transformed Velocity')
        axs[1].axis('equal')
        axs[1].set_xlim(xl,xu)
        axs[1].set_ylim(yl,yu)
        plt.tight_layout()

def plot_pressure(Z, Cp, circle, Z_transformed, Cp_transformed, airfoil):

    # Set the plot limits --- Change these if needed
    yl = -3
    yu = 3
    xl = -3
    xu = 3

    fig = plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.05], hspace=0.3)

    # Create subplots for the pressure plots
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    cp_levels = np.linspace(-1.0, 1.0, 100)

    # Plot the original pressure
    cp_contour_0 = axs[0].contourf(Z.real, Z.imag, Cp, levels=cp_levels, extend='both')
    axs[0].plot(np.real(circle), np.imag(circle), color='red')
    axs[0].set_title('Original Pressure')
    axs[0].axis('equal')
    axs[0].set_xlim(xl, xu)
    axs[0].set_ylim(yl, yu)

    # Plot the transformed pressure
    cp_contour_1 = axs[1].contourf(Z_transformed.real, Z_transformed.imag, Cp_transformed, levels=cp_levels, extend='both')
    axs[1].plot(np.real(airfoil), np.imag(airfoil), color='red')
    axs[1].set_title('Transformed Pressure')
    axs[1].axis('equal')
    axs[1].set_xlim(xl, xu)
    axs[1].set_ylim(yl, yu)

    # Add a colorbar below the plots
    cbar_ax = fig.add_subplot(gs[1, :])
    cbar = fig.colorbar(cp_contour_0, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Cp')
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])

def plot_pressure_distribution(X, Cp, X_transformed, Cp_transformed):
    # Same as plot_pressure_difference but without the lift calculation
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    axs[0].scatter(X, Cp, s=8)
    axs[0].set_title('Original Pressure Distribution')
    axs[1].scatter(X_transformed, Cp_transformed, s=8)
    axs[1].set_title('Transformed Pressure Distribution')
    plt.tight_layout()

def plot_grid(Z, circle, Z_transformed, airfoil):
    


    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    axs[0].scatter(Z.real, Z.imag, s=1)
    axs[0].plot(np.real(circle), np.imag(circle), color='red')
    axs[0].set_title('Original Grid')

    axs[1].scatter(Z_transformed.real, Z_transformed.imag, s=1)
    axs[1].plot(np.real(airfoil), np.imag(airfoil), color='red')
    axs[1].set_title('Transformed Grid')
    axs[1].axis('equal')

    plt.tight_layout()

def plot_pressure_difference(Z,Cp, Z_transformed, Cp_transformed,yc):    
    Xc_upper, Xc_lower, Cpc_upper, Cpc_lower = airfoil_split(Z,Cp,Z,yc)
    L_circle = calculate_lift(Z, Cp,Z,yc)

    Xa_upper, Xa_lower, Cpa_upper, Cpa_lower = airfoil_split(Z_transformed,Cp_transformed,Z,yc)
   
    # Thresholding to remove singularities
    threshold = 1e6
    Xa_lower, Cpa_lower = rem_threshold(Cpa_lower, Xa_lower, threshold)
    Xa_upper, Cpa_upper = rem_threshold(Cpa_upper, Xa_upper, threshold)

    L_airfoil = calculate_lift(Z_transformed, Cp_transformed,Z,yc)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    
    axs[0].scatter(Xc_upper, Cpc_upper, s=8, label='Upper Surface', color='red')
    axs[0].scatter(Xc_lower, Cpc_lower, s=8, label='Lower Surface', color='blue')
    axs[0].set_title('Lift (Circle) = {:.4f}'.format(L_circle))
    axs[0].set_xlabel('X')
    axs[0].set_ylabel('Cp')
    axs[1].scatter(Xa_upper, Cpa_upper, s=8, label='Upper Surface', color='red')
    axs[1].scatter(Xa_lower, Cpa_lower, s=8, label='Lower Surface', color='blue')
    axs[1].set_title('Lift (Airfoil) = {:.4f}'.format(L_airfoil))
    axs[1].set_xlabel('X')
    axs[1].set_ylabel('Cp')
    plt.tight_layout()

def rem_threshold(Cp,X, threshold):
    rem_indices = np.where(abs(Cp) > threshold)
    X = np.delete(X, rem_indices)
    Cp = np.delete(Cp, rem_indices)
    return X, Cp

def plot_fitting(airfoil_coords, airfoil_fit):
    plt.figure(figsize=(8, 6))
    plt.plot(np.real(airfoil_coords), np.imag(airfoil_coords), 'ro', label='Airfoil Data')
    plt.plot(np.real(airfoil_fit), np.imag(airfoil_fit), 'b-', label='Fitted Airfoil')
    plt.axis('equal')
    plt.legend()
    plt.title('Joukowski Transform Fitting')
    plt.xlabel('Re(ζ)')
    plt.ylabel('Im(ζ)')
    plt.grid()

def plot_all(Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, yc, stride=5, AoA=0):
    fig, axs = plt.subplots(3, 2, figsize=(8, 9))
    plt.suptitle('Joukowski Transformation Results', fontsize=16)

    # Grid Plot
    #axs[0, 0].scatter(Z.real, Z.imag, s=1)
    axs[0, 0].plot(np.real(circle), np.imag(circle), color='red')
    axs[0, 0].set_title('Original Grid')
    axs[0, 0].axis('equal')

    #axs[0, 1].scatter(Z_transformed.real, Z_transformed.imag, s=1)
    axs[0, 1].plot(np.real(airfoil), np.imag(airfoil), color='red')
    axs[0, 1].set_title('Transformed Grid')
    axs[0, 1].axis('equal')

    # Potential Plot
    axs[1, 0].contour(Z.real, Z.imag, Phi, levels=50)
    axs[1, 0].plot(np.real(circle), np.imag(circle), color='red')
    axs[1, 0].set_title('Original Velocity Field')
    axs[1, 0].axis('equal')

    axs[1, 1].contour(Z_transformed.real, Z_transformed.imag, Phi, levels=50)
    axs[1, 1].plot(np.real(airfoil), np.imag(airfoil), color='red')
    axs[1, 1].set_title('Transformed Velocity Field')
    axs[1, 1].axis('equal')

    # Pressure Difference Plot
    Xc_upper, Xc_lower, Cpc_upper, Cpc_lower = airfoil_split(Z,Cp,Z,yc)
    L_circle = calculate_lift(Z, Cp,Z,yc)

    Xa_upper, Xa_lower, Cpa_upper, Cpa_lower = airfoil_split(Z_transformed,Cp_transformed,Z,yc)
    L_airfoil = calculate_lift(Z_transformed, Cp_transformed,Z,yc)

    #Xa_lower = Xa_lower[:,1:-2]
    #Xa_upper = Xa_upper[:,1:-2]
    #Cpa_upper = Cpa_upper[:,1:-2]
    #Cpa_lower = Cpa_lower[:,1:-2]

    # Thresholding to remove singularities
    threshold = 1e6
    Xa_lower, Cpa_lower = rem_threshold(Cpa_lower, Xa_lower, threshold)
    Xa_upper, Cpa_upper = rem_threshold(Cpa_upper, Xa_upper, threshold)


    axs[2, 0].scatter(Xc_upper, Cpc_upper, s=8, label='Upper Surface', color='red')
    axs[2, 0].scatter(Xc_lower, Cpc_lower, s=8, label='Lower Surface', color='blue')
    axs[2, 0].set_title('Lift (Circle) = {:.4f}'.format(L_circle))

    axs[2, 1].scatter(Xa_upper, Cpa_upper, s=8, label='Upper Surface', color='red')
    axs[2, 1].scatter(Xa_lower, Cpa_lower, s=8, label='Lower Surface', color='blue')
    axs[2, 1].set_title('Lift (Airfoil) = {:.4f}'.format(L_airfoil))

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit suptitle


def JT_solver(x_in, y_in, R_in, AoA_in, U_in):
    # Circle parameters
    xc = x_in
    yc = y_in
    a =  R_in

    # Grid parameters
    N_r = 100
    N_theta = 145
    
    # Joukowski Transformation parameters
    c = 1
    # Flow parameters
    AoA = AoA_in*np.pi/180
    U_inf = U_in

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

    return Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, L_circle, L_airfoil, AoA, u_z, v_z, u_transformed, v_transformed


def interactive_plot():
    def plot_all_dynamic(Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, yc, pressure_flag, velocity_flag, u_z, v_z, u_transformed, v_transformed,xc):
        # Clear previous plots
        axs[0, 0].cla()
        axs[0, 1].cla()
        axs[1, 0].cla()
        axs[1, 1].cla()
        axs[2, 0].cla()
        axs[2, 1].cla()

        
        yl = -3#yc-3
        yu = 3#yc+3
        xl = -3#xc-3
        xu = 3#xc+3

        # Plot Grids
        axs[0, 0].plot(np.real(circle), np.imag(circle), color='red')
        axs[0, 0].set_title('Original Grid')
        axs[0, 0].axis('equal')
        axs[0, 0].set_xlim(xl, xu)
        axs[0, 0].set_ylim(yl, yu)

        axs[0, 1].plot(np.real(airfoil), np.imag(airfoil), color='red')
        axs[0, 1].set_title('Transformed Grid')
        axs[0, 1].axis('equal')
        axs[0, 1].set_xlim(xl, xu)
        axs[0, 1].set_ylim(yl, yu)

        # Conditional plotting based on flags
        if pressure_flag == "Pressure Dist.":
            
            # Pressure Difference Plot
            Xc_upper, Xc_lower, Cpc_upper, Cpc_lower = airfoil_split(Z, Cp, Z, yc)
            L_circle = calculate_lift(Z, Cp, Z, yc)

            Xa_upper, Xa_lower, Cpa_upper, Cpa_lower = airfoil_split(Z_transformed, Cp_transformed, Z, yc)
            L_airfoil = calculate_lift(Z_transformed, Cp_transformed, Z, yc)

            axs[2,0].cla()
            axs[2, 0].scatter(Xc_upper, Cpc_upper, s=8, label='Upper Surface', color='red')
            axs[2, 0].scatter(Xc_lower, Cpc_lower, s=8, label='Lower Surface', color='blue')
            axs[2, 0].set_title(f'Lift (Circle) = {L_circle:.4f}')
            axs[2, 0].set_xlim(xl, xu)

            axs[2,1].cla()
            axs[2, 1].scatter(Xa_upper, Cpa_upper, s=8, label='Upper Surface', color='red')
            axs[2, 1].scatter(Xa_lower, Cpa_lower, s=8, label='Lower Surface', color='blue')
            axs[2, 1].set_title(f'Lift (Airfoil) = {L_airfoil:.4f}')
            axs[2, 1].set_xlim(xl, xu)
        
        elif pressure_flag == "Pressure":
            # Set the plot limits
            cp_levels = np.linspace(-1.0, 1.0, 100)

            # Pressure Plots
            axs[2, 0].cla()
            cp_contour_0 = axs[2, 0].contourf(Z.real, Z.imag, Cp, levels=cp_levels, extend='both')
            axs[2, 0].plot(np.real(circle), np.imag(circle), color='red')
            axs[2, 0].set_title('Original Pressure')
            axs[2, 0].axis('equal')
            axs[2, 0].set_xlim(xl, xu)
            axs[2, 0].set_ylim(yl, yu)

            axs[2, 1].cla()
            cp_contour_1 = axs[2, 1].contourf(Z_transformed.real, Z_transformed.imag, Cp_transformed, levels=cp_levels, extend='both')
            axs[2, 1].plot(np.real(airfoil), np.imag(airfoil), color='red')
            axs[2, 1].set_title('Transformed Pressure')
            axs[2, 1].axis('equal')
            axs[2, 1].set_xlim(xl, xu)
            axs[2, 1].set_ylim(yl, yu)

        if velocity_flag == "Velocity":
            # Velocity Plots
            stride = 10  # You can change this for better resolution
            s1 = 'default'
            s2 = 20

            axs[1, 0].cla()
            axs[1, 0].quiver(Z.real[::stride, ::stride], Z.imag[::stride, ::stride], u_z[::stride, ::stride], v_z[::stride, ::stride])
            axs[1, 0].plot(np.real(circle), np.imag(circle), color='red')
            axs[1, 0].set_title('Original Velocity')
            axs[1, 0].axis('equal')
            axs[1, 0].set_xlim(xl, xu)
            axs[1, 0].set_ylim(yl, yu)

            axs[1,1].cla()
            axs[1, 1].quiver(Z_transformed.real[::stride, ::stride], Z_transformed.imag[::stride, ::stride], u_transformed[::stride, ::stride], v_transformed[::stride, ::stride], scale=s2)
            axs[1, 1].plot(np.real(airfoil), np.imag(airfoil), color='red')
            axs[1, 1].set_title('Transformed Velocity')
            axs[1, 1].axis('equal')
            axs[1, 1].set_xlim(xl, xu)
            axs[1, 1].set_ylim(yl, yu)
        
        elif velocity_flag == "Potential":
            # Potential Plot
            axs[1, 0].cla()  # Clear only the current contour plot
            axs[1, 0].contour(Z.real, Z.imag, Phi, levels=50)
            axs[1, 0].plot(np.real(circle), np.imag(circle), color='red')
            axs[1, 0].set_title('Original Potential')
            axs[1, 0].axis('equal')
            axs[1, 0].set_xlim(xl, xu)
            axs[1, 0].set_ylim(yl, yu)

            axs[1, 1].cla()  # Clear only the current contour plot
            axs[1, 1].contour(Z_transformed.real, Z_transformed.imag, Phi, levels=50)
            axs[1, 1].plot(np.real(airfoil), np.imag(airfoil), color='red')
            axs[1, 1].set_title('Transformed Potential')
            axs[1, 1].axis('equal')
            axs[1, 1].set_xlim(xl, xu)
            axs[1, 1].set_ylim(yl, yu)

        plt.tight_layout(rect=[0, 0, 0.75, 1])
        return axs

    def update(val):
        # Get slider values
        x_in = sx.val
        y_in = sy.val
        R_in = sr.val
        AoA_in = sAoA.val
        U_in = sU.val
        pressure_flag = pressure_selector.value_selected
        velocity_flag = velocity_selector.value_selected

        # Solve for new values
        Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, _, _, _, u_z, v_z, u_transformed, v_transformed = JT_solver(x_in, y_in, R_in, AoA_in, U_in)

        # Update the plots without clearing the axes
        plot_all_dynamic(Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, y_in, pressure_flag, velocity_flag, u_z, v_z, u_transformed, v_transformed,x_in)

        # Update the text boxes with the new slider values
        txt_x.set_val(f"{x_in:.2f}")
        txt_y.set_val(f"{y_in:.2f}")
        txt_r.set_val(f"{R_in:.2f}")
        txt_AoA.set_val(f"{AoA_in:.2f}")
        txt_U.set_val(f"{U_in:.2f}")
        
        fig.canvas.draw_idle()

    def update_text(val, slider):
        try:
            new_val = float(val)
            slider.set_val(new_val)
        except ValueError:
            pass  # Ignore invalid input

    # Initial Input
    x_in = 0
    y_in = 0
    R_in = 1
    AoA_in = 20
    U_in = 1
    pressure_flag = "Pressure"
    velocity_flag = "Potential"  # Set initial flag to "Velocity"

    Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, _, _, _, u_z, v_z, u_transformed, v_transformed = JT_solver(x_in, y_in, R_in, AoA_in, U_in)
    fig, axs = plt.subplots(3, 2, figsize=(10, 9))
    plt.suptitle('Joukowski Transformation Results', fontsize=16)

    # Create the initial plot
    plot_all_dynamic(Z, circle, Z_transformed, airfoil, Phi, Cp, Cp_transformed, y_in, pressure_flag, velocity_flag, u_z, v_z, u_transformed, v_transformed,x_in)

    # Adjust the layout for the sliders
    plt.subplots_adjust(right=0.75)  # Increase the space for the sliders (0.85)

    # Slider Axes (Positioned on the right)
    x_pos = 0.8
    y_pos = 0.4
    axcolor = 'lightgoldenrodyellow'
    ax_x = plt.axes([x_pos, y_pos, 0.12, 0.02], facecolor=axcolor)
    ax_y = plt.axes([x_pos, y_pos+0.03, 0.12, 0.02], facecolor=axcolor)
    ax_r = plt.axes([x_pos, y_pos+0.06, 0.12, 0.02], facecolor=axcolor)
    ax_AoA = plt.axes([x_pos, y_pos+0.09, 0.12, 0.02], facecolor=axcolor)
    ax_U = plt.axes([x_pos, y_pos+0.12, 0.12, 0.02], facecolor=axcolor)

    # Create Text Boxes for each input
    txt_x = TextBox(plt.axes([x_pos+0.04, y_pos+0.15, 0.05, 0.02]), 'x ', initial=str(f"{x_in:.2f}"))
    txt_y = TextBox(plt.axes([x_pos+0.04, y_pos+0.18, 0.05, 0.02]), 'y ' , initial=str(f"{y_in:.2f}"))
    txt_r = TextBox(plt.axes([x_pos+0.04, y_pos+0.21, 0.05, 0.02]), 'R ', initial=str(f"{R_in:.2f}"))
    txt_AoA = TextBox(plt.axes([x_pos+0.04, y_pos+0.24, 0.05, 0.02]), 'AoA ', initial=str(f"{AoA_in:.2f}"))
    txt_U = TextBox(plt.axes([x_pos+0.04, y_pos+0.27, 0.05, 0.02]), 'U ', initial=str(f"{U_in:.2f}"))

    # Link text boxes to sliders
    txt_x.on_text_change(lambda text: update_text(text, sx))
    txt_y.on_text_change(lambda text: update_text(text, sy))
    txt_r.on_text_change(lambda text: update_text(text, sr))
    txt_AoA.on_text_change(lambda text: update_text(text, sAoA))
    txt_U.on_text_change(lambda text: update_text(text, sU))

    # Create Sliders
    sx = Slider(ax_x, 'x', -2.0, 2.0, valinit=x_in)
    sy = Slider(ax_y, 'y ', -2.0, 2.0, valinit=y_in)
    sr = Slider(ax_r, 'R', 0.5, 2.0, valinit=R_in)
    sAoA = Slider(ax_AoA, 'AoA', 0.0, 45.0, valinit=AoA_in)
    sU = Slider(ax_U, 'U', 0.1, 5.0, valinit=U_in)

    # Create Flag Selectors for Pressure and Velocity using Radio Buttons
    pressure_selector = RadioButtons(plt.axes([x_pos, y_pos - 0.25 , 0.15, 0.05]), ['Pressure', 'Pressure Dist.'])
    velocity_selector = RadioButtons(plt.axes([x_pos, y_pos - 0.15, 0.15, 0.05]), ['Potential', 'Velocity'])

    # Update Function for Sliders
    sx.on_changed(update)
    sy.on_changed(update)
    sr.on_changed(update)
    sAoA.on_changed(update)
    sU.on_changed(update)
    pressure_selector.on_clicked(update)
    velocity_selector.on_clicked(update)

    plt.show()