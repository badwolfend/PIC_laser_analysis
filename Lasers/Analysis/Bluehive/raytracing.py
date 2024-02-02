import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def plot_all_rays(ax, paths):

    # Density plot
    # Plot the density distribution
    plt.imshow(Z, extent=[min(x), max(x), min(y), max(y)], origin='lower', cmap='Greys', alpha=1, aspect='auto')

    # Ray's path
    for path in paths:
        plt.plot(path[:, 0], path[:, 1], 'xkcd:green', linewidth=4, label='Ray')

        # Adding an arrow at the end of the ray's path
        end_point = path[-1]
        direction = path[-1] - path[-2]  # Direction vector for the arrow
        arrow_length = np.linalg.norm(direction)  # Length of the arrow
        arrow_length_scaled = arrow_length * 5  # Scale up for visibility

        # Ensure the arrow is visible and properly scaled
        if arrow_length_scaled > 0:
            plt.arrow(path[-2, 0], path[-2, 1], direction[0], direction[1], 
                    head_width=0.4, head_length=0.4, fc='xkcd:green', ec='xkcd:green')
        
    # Limite plot to the region of interest
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plasma Density and Ray Path')
    plt.legend()
    ax.set_aspect('equal')
    plt.show()

def plot_ray_path(circle=False):
    # Plotting
    fig = plt.figure(figsize=(8, 6))

    # Set aspect ratio of plot to be auto
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')

    # Density plot
    # Plot the density distribution
    plt.imshow(Z, extent=[min(x), max(x), min(y), max(y)], origin='lower', cmap='Greys', alpha=1, aspect='auto')

    if circle:
        # Find the location where path is closest to the center
        min_distance = np.inf
        min_index = 0
        for i, point in enumerate(path):
            distance = np.linalg.norm(point)
            if distance < min_distance:
                min_distance = distance
                min_index = i
        print(f"Closest approach to center: {min_distance} at point {path[min_index]}")
        closest = path[min_index]

        # Define the circle
        # Specify the center (x, y) and the radius
        center = (0, 0)  # Example center point
        radius = np.sqrt(closest[0]**2+closest[1]**2)  # Example radius
        neturn = plasma_density(np.array(0), np.array(radius), A, B)
        print(f"Computed impact paramter (b/L): {radius*np.sqrt(1-neturn)}")

        # Create a circle patch with dotted edge
        circle = Circle(center, radius, fill=False, linestyle='--', linewidth=4, edgecolor='xkcd:salmon', label='Closest approach to center')

        # Add the circle to the plot
        ax.add_patch(circle)
        # contour = plt.contourf(X, Y, Z, levels=100, cmap='Greys')
        # plt.colorbar(contour)

    # Ray's path
    plt.plot(path[:, 0], path[:, 1], 'xkcd:green', linewidth=4, label='Ray')

    # Adding an arrow at the end of the ray's path
    end_point = path[-1]
    direction = path[-1] - path[-2]  # Direction vector for the arrow
    arrow_length = np.linalg.norm(direction)  # Length of the arrow
    arrow_length_scaled = arrow_length * 5  # Scale up for visibility

    # Ensure the arrow is visible and properly scaled
    if arrow_length_scaled > 0:
        plt.arrow(path[-2, 0], path[-2, 1], direction[0], direction[1], 
                head_width=0.4, head_length=0.4, fc='xkcd:green', ec='xkcd:green')
        
    # Limite plot to the region of interest
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plasma Density and Ray Path')
    plt.legend()
    ax.set_aspect('equal')
    plt.show()

def plasma_density(x, y, A=1, B=0.01, C=0, type='exponential'):
    """
    Computes the plasma density at a given point (x, y).
    A and B are parameters that define the density distribution.
    """

    if type == 'linear':
        if (x.size) <= 1:
            # Compute the linear density distribution ramping from 0 to A for x < B
            dens = A * ((x-C)/(B-C))
            if x>=(B):
                dens = A
            if x<=C:
                dens = 0
        else:
            # Compute the 2D linear density distribution ramping from 0 to A for x < B and uniform in y
            dens = np.zeros_like(x) 
            dens = A * ((x-C)/(B-C))
            dens[x>=(B)] = A
            dens[x<=C] = 0

    if type == 'exponential':
        # Compute the exponential density distribution for radius >= B
        if (x.size) <= 1:
            r = x**2+y**2
            if r>B:
                dens = A * np.exp(-1 * (r-B))
            else:
                dens = A 
        else:
            r2 = x**2 + y**2
            r = np.sqrt(r2)
            dens = np.ones_like(r2)
            dens[r>=B] = A * np.exp(-1 * (r[r>=B]-B))
        
    return dens

def compute_gradient(x, y, A=1, B=0.01, C=0, type='exponential', h=1e-5):
    """
    Computes the gradient of the plasma density at (x, y) using finite differences.
    h is the step size for the finite difference calculation.
    """
    density_x1 = plasma_density(x + h, y, A, B, C, type)
    density_x0 = plasma_density(x - h, y, A, B, C, type)
    density_y1 = plasma_density(x, y + h, A, B, C, type)
    density_y0 = plasma_density(x, y - h, A, B, C, type)
    
    grad_x = (density_x1 - density_x0) / (2 * h)
    grad_y = (density_y1 - density_y0) / (2 * h)
    
    return -np.array([grad_x, grad_y])

def rk4_step(f, r, k, dt, A=1, B=0.01, C=0, type='exponential'):
    """
    Performs a single step of the RK4 integration for 2D vectors, taking into account plasma density.
    f: Function to integrate, returning the gradient vector.
    r, k: Current position and wave vector (direction) vectors in 2D.
    dt: Step size.
    A, B: Parameters for the plasma density distribution.
    """
    c2oomega = (clight**2/V0**2)/omega

    k1_r =  c2oomega*k*dt
    k1_k =  0.5*omega*f(r[0], r[1], A, B, C, type) * dt
    
    k2_r = c2oomega*(k + 0.5*k1_k) * dt
    k2_k = 0.5*omega*f(r[0] + 0.5*k1_r[0], r[1] + 0.5*k1_r[1], A, B, C, type) * dt
    
    k3_r = c2oomega*(k + 0.5*k2_k) * dt
    k3_k = 0.5*omega*f(r[0] + 0.5*k2_r[0], r[1] + 0.5*k2_r[1], A, B, C, type) * dt
    
    k4_r = c2oomega*(k + k3_k) * dt
    k4_k = 0.5*omega*f(r[0] + k3_r[0], r[1] + k3_r[1], A, B, C, type) * dt
    
    r_next = r + (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6
    v_next = k + (k1_k + 2*k2_k + 2*k3_k + k4_k) / 6
    
    return r_next, v_next

def ray_trace(r, k, dt, steps, A=1, B=0.01, C=0, type='exponential'):
    """
    Ray tracing in a 2D plasma density gradient with specified density distribution.
    r: Initial position of the ray (2D vector).
    k: Initial wave vector (direction) of the ray (2D vector).
    ds: Step size for integration.
    steps: Number of integration steps.
    A, B: Parameters for the plasma density distribution.
    """
    path = [r]
    
    for _ in range(steps):
        r, k = rk4_step(compute_gradient, r, k, dt, A, B, C, type)
        path.append(r)
        
    return np.array(path)


def test_exponential_plasma():

    # Generate a meshgrid for density plot
    x = np.linspace(-4, 4, 1000)
    y = np.linspace(0, 4, 1000)
    X, Y = np.meshgrid(x, y)
    A = 1.0
    B = 0.4
    C = 0
    type = 'exponential'
    Z = plasma_density(X, Y, A=A, B=B)  # Plasma density function

    # Laser parameters
    clight = 3e8  # Speed of light in m/s
    lambda_laser = 532e-9  # Wavelength in meters
    omega_laser = 2*np.pi*clight/lambda_laser  # Angular frequency (rad/s)

    # Define the nondimensionalization 
    L0 = 10*lambda_laser  # Length scale in meters
    V0 = clight  # Time scale in seconds
    T0 = L0/V0  # Velocity scale in m/s

    # Now none-dimensionalize the parameters
    omega = omega_laser*T0  # Nondimensionalized angular frequency
    r0 = np.array([-4, 0.5])  # Starting point of the ray (x, y) already nondimensionalized
    direction = np.array([1, 0])  # Initial direction
    k0 = L0*direction/np.linalg.norm(direction)*omega_laser/clight  # Initial nondimensionalized direction

    # Compute the ray's path
    dt = 0.01  # Step size
    steps = 600  # Number of steps

    return A, B, C, type, x, X, y, Y, Z, dt, steps, L0, V0, T0, omega, clight, lambda_laser, r0, k0


# Assuming previous definitions for plasma_density, compute_gradient, rk4_step, and ray_trace are available
# Laser parameters
clight = 3e8  # Speed of light in m/s
lambda_laser = 532e-9  # Wavelength in meters
omega_laser = 2*np.pi*clight/lambda_laser  # Angular frequency (rad/s)

# Define the nondimensionalization 
L0 = 1*lambda_laser  # Length scale in meters
V0 = clight  # Time scale in seconds
T0 = L0/V0  # Velocity scale in m/s

# Generate a meshgrid for density plot
x = np.linspace(0, 12, 1000)
y = np.linspace(0, 12, 1000)
X, Y = np.meshgrid(x, y)
A = 1.0
C = 1.0
B = C+1.5*lambda_laser/L0
type = 'linear'

Z = plasma_density(X, Y, A=A, B=B, C=C, type=type)  # Plasma density function

# Now none-dimensionalize the parameters
omega = omega_laser*T0  # Nondimensionalized angular frequency
r0 = np.array([0, 0])  # Starting point of the ray (x, y) already nondimensionalized
direction = np.array([1, 1])  # Initial direction
k0 = L0*direction/np.linalg.norm(direction)*omega_laser/clight  # Initial nondimensionalized direction

# Compute the ray's path
dt = 0.01  # Step size
steps = 600  # Number of steps
if False:
    A, B, C, type, x, X, y, Y, Z, dt, steps, L0, V0, T0, omega, clight, lambda_laser, r0, k0,  = test_exponential_plasma()
    path = ray_trace(r0, k0, dt, steps, A, B, C, type)

    # Now plot things 
    plot_ray_path(circle=True)

else:
    Nr = 10
    DR =  np.max(y)/2
    # Plotting
    fig = plt.figure(figsize=(8, 6))

    # Set aspect ratio of plot to be auto
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    paths = []
    for i in range(Nr):
        r0 = np.array([0, DR*(i)/Nr])
        paths.append(ray_trace(r0, k0, dt, steps, A, B, C, type))
    plot_all_rays(ax, paths)
