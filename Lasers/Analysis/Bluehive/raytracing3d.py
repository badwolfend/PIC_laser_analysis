import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

def find_critical_dens(wl):
    """
    - wl in micron #
    - dens in /m^3   
    """
 
    n_crit = (10**6)*(1.1e21)/(wl**2)
    return n_crit

def plot_ray_path_3d_pyvis(plotter, path, plasma_density_3d=None, A=1, B=0.01, C=0, type='bennett', x=None, y=None, z=None):
    N=1
    X, Y, Z = np.meshgrid(x[::N], y[::N], z[::N], indexing='ij')  # This Z is different from plasma density Z, consider renaming for clarity

    # Compute the density at each point in the meshgrid
    dens = plasma_density_3d(X, Y, Z, A, B, C, type)+0.001

    # Convert the numpy arrays into a PyVista grid
    grid = pv.StructuredGrid(X, Y, Z)
    grid["densities"] = dens.flatten(order="F") # Add the density data to the grid

    # Plot isosurfaces
    if plotter == None:
      plotter = pv.Plotter(off_screen=True)
      # plotter = pv.Plotter()
    if True:
      threshed_mesh = grid.threshold(value=0.01, scalars="densities")
      slices = threshed_mesh.slice_orthogonal()
      # plotter.add_mesh(threshed_mesh, opacity="geom_r", cmap='reds', clim=(0, 1))  # Change 'cyan' to your desired color
      plotter.add_mesh(slices, opacity=1.0, cmap='terrain', style="points", point_size=20, render_points_as_spheres=True, interpolate_before_map=True)  # Change 'cyan' to your desired color
      # plotter.set_viewup([1,0,0])

      # plotter.set_focus(threshed_mesh.points[1000])

      # plotter.add_mesh(slices, opacity=1.0, cmap='terrain',  interpolate_before_map=True, smooth_shading=True)  # Change 'cyan' to your desired color

      # plotter.add_mesh_threshold(grid, scalars="densities", invert=True, cmap="reds", opacity=0.5)
    
    else:
      # Calculate the azimuthal angle for each point in the grid
      points = grid.points
      x, y, z = points[:, 0], points[:, 1], points[:, 2]
      theta = np.arctan2(y, x)  # Arctan2 returns angles in the range [-pi, pi]

      # Define your angle range in radians (e.g., -pi/4 to pi/4)
      angle_min = -np.pi+3*np.pi/4
      # angle_min = -np.pi
      angle_max = np.pi

      # Create a mask for points within the specified angle range
      angle_mask = (theta >= angle_min) & (theta <= angle_max)

      # Apply the mask to the densities (set densities outside the range to NaN or another value)
      filtered_densities = np.copy(grid["densities"])
      filtered_densities[~angle_mask] = np.nan  # Using NaN to hide points outside the angle range

      # Update the grid with the filtered densities
      grid['filtered_densities'] = filtered_densities
      threshed_mesh = grid.threshold(value=0.01, scalars="filtered_densities")
      # plotter.add_mesh(threshed_mesh, opacity=1, cmap='terrain', style="points", point_size=20, render_points_as_spheres=True, interpolate_before_map=True, smooth_shading=True, lighting=True, metallic=0.7, roughness=0.7 )  # Change 'cyan' to your desired color
      plotter.add_mesh(threshed_mesh, cmap='terrain', specular=1.0, interpolate_before_map=True, smooth_shading=True, lighting=True, metallic=0.7, roughness=0.7 )  # Change 'cyan' to your desired color

      # plotter.add_mesh_threshold(grid, scalars="filtered_densities", invert=False, cmap="reds")
      # Optionally, set the scalar bar to ignore NaN values
      plotter.scalar_bar.SetLookupTable(plotter.mapper.GetLookupTable())
      plotter.scalar_bar.GetLookupTable().SetNanColor(0, 0, 0, 0)  # Set NaN color to transparent (RGBA)


    # Add the ray path to the plotter
    plotter.add_mesh(path, color='green', line_width=2, label='Ray Path')

    # Calculate the direction vector for the arrow
    direction_vector = path[-1] - path[-2]
    direction_vector = direction_vector / np.linalg.norm(direction_vector)  # Normalize the vector

    # The position where the arrow will be placed (the end of the ray path)
    arrow_position = path[-1]

    # Parameters for the arrow
    scale = 0.3  # Adjust the scale of the arrow as needed
    arrow_start = arrow_position - direction_vector * scale  # Adjust start position based on scale
    arrow_start = arrow_position

    # Create the arrow
    arrow = pv.Arrow(start=arrow_start, direction=direction_vector, scale=scale)

    # Add the arrow to the plotter
    plotter.add_mesh(arrow, color='red', line_width=4, label='Arrow')
    # plotter.show()
    return plotter


def plot_ray_path_3d(path, plasma_density_3d=None, A=1, B=0.01, C=0, type='bennett', X=None, Y=None, Z=None):
    """
    Plots the ray path in a 3D plasma density distribution.
    
    Args:
    - path: Array of points representing the ray path in 3D.
    - plasma_density_3d: A function to calculate plasma density (optional for visualization).
    - A, B: Parameters for the plasma density function.
    - x, y, z: Meshgrids for plasma density (if plasma_density_3d is provided).
    - circle: Not directly applicable in 3D as in 2D, but could represent a sphere for the closest approach.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # If a density function and meshgrids are provided, attempt a simple density visualization
    if plasma_density_3d and X is not None and Y is not None and Z is not None:
        if True:
          N=3
          X, Y, Z = np.meshgrid(x[::N], y[::N], z[::N])  # This Z is different from plasma density Z, consider renaming for clarity

          # Compute the density at each point in the meshgrid
          dens = plasma_density_3d(X, Y, Z, A, B, C, type)
          xflat = X.flatten()
          yflat = Y.flatten()
          zflat = Z.flatten()
          densflat = dens.flatten()
          scatter = ax.scatter(xflat, yflat, zflat, c=densflat, cmap='Reds', alpha=0.5, s=1)
          # scatter = ax.scatter(X.flatten()[::N], Y.flatten()[::N], Z.flatten()[::N], c=dens.flatten()[::N], cmap='Reds', alpha=0.1)
          fig.colorbar(scatter, ax=ax, label='Density')
        else:
          N=3
          X, Y, Z = np.meshgrid(x[::N], y[::N], z[::N])  # This Z is different from plasma density Z, consider renaming for clarity

          # Compute the density at each point in the meshgrid
          dens = plasma_density_3d(X, Y, Z, A, B)
        
          # Convert the numpy arrays into a PyVista grid
          grid = pv.StructuredGrid(X, Y, Z)
          grid["densities"] = densities.flatten(order='F')  # Add the density data to the grid

          # Plot isosurfaces
          plotter = pv.Plotter()
          plotter.add_mesh_threshold(grid, scalars="densities", invert=True)
          plotter.show()
    else:
        # This is a placeholder: real 3D density visualization may require volume rendering or isosurfaces
        # For simplicity, you might just plot a subset of points or use contours at certain slices
        pass

    # Plotting the ray path
    ax.plot3D(path[:, 0], path[:, 1], path[:, 2], color='green', linewidth=8, label='Ray Path')


    # Adding an arrow at the end of the ray's path in 3D
    end_point = path[-1]
    direction = path[-1] - path[-2]  # Direction vector for the arrow

    # Normalize the direction for the arrow
    direction_norm = direction / np.linalg.norm(direction)

    # Scaling factor for the arrow length, adjust as necessary for visibility
    arrow_length_scaled = 0.2 # Adjust based on your coordinate system scale

    # Calculate the arrow's end point (start + scaled direction)
    arrow_end = end_point + direction_norm * arrow_length_scaled

    # Use quiver to add the arrow
    # Note: quiver in 3D doesn't support arrowheads, so this will just draw a line in the direction
    ax.quiver(end_point[0], end_point[1], end_point[2], 
              direction_norm[0], direction_norm[1], direction_norm[2], 
              length=arrow_length_scaled, color='xkcd:green', linewidth=8)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D Plasma Density and Ray Path')
    plt.legend()
    plt.show()

def plasma_density_3d(x, y, z, A=1, B=0.01, C=0, type='exponential'):
    """
    Computes the plasma density at a given point (x, y, z) in 3D space.
    A and B are parameters that define the density distribution.
    """
  

    if type == 'exponential':
        # Compute the exponential density distribution for radius >= B in 3D
        if np.isscalar(x):
            r = np.sqrt(x**2 + y**2 + z**2)
            dens = A * np.exp(-1 * (r-B)) if r > B else A
        else:
            r2 = x**2 + y**2 + z**2
            r = np.sqrt(r2)
            dens = np.ones_like(r2) * A
            dens[r > B] = A * np.exp(-1 * (r[r > B] - B))

    if type == 'bennett':
        # Compute the Bennett pinch density distribution in 3D
        if np.isscalar(x):
            rc = np.sqrt(x**2 + y**2)
            dens = A*(1-(rc**2 / B**2)) if rc <= B else 0
        else:
            r2 = x**2 + y**2
            r = np.sqrt(r2)
            dens = A*(1-(r2 / B**2))
            dens[r > B] = 0
            
    return dens

def compute_gradient_3d(x, y, z, A=1, B=0.01, C=0, type='exponential', h=1e-5):
    """
    Computes the gradient of the plasma density at (x, y, z) using finite differences in 3D.
    h is the step size for the finite difference calculation.
    """
    # Compute density for small increments in each direction
    density_x1 = plasma_density_3d(x + h, y, z, A, B, C, type)
    density_x0 = plasma_density_3d(x - h, y, z, A, B, C, type)
    density_y1 = plasma_density_3d(x, y + h, z, A, B, C, type)
    density_y0 = plasma_density_3d(x, y - h, z, A, B, C, type)
    density_z1 = plasma_density_3d(x, y, z + h, A, B, C, type)
    density_z0 = plasma_density_3d(x, y, z - h, A, B, C, type)
    
    # Calculate gradient components
    grad_x = (density_x1 - density_x0) / (2 * h)
    grad_y = (density_y1 - density_y0) / (2 * h)
    grad_z = (density_z1 - density_z0) / (2 * h)
    
    return -np.array([grad_x, grad_y, grad_z])

def rk4_step_3d(f, r, k, dt, A=1, B=0.01, C=0, type='exponential'):
    """
    Performs a single step of the RK4 integration for 3D vectors, taking into account plasma density.
    f: Function to integrate, returning the gradient vector.
    r, k: Current position and wave vector (direction) vectors in 3D.
    dt: Step size.
    A, B, C: Parameters for the plasma density distribution.
    """
    # Calculate the conversion factor for wave vector based on medium properties
    c2oomega = (clight**2 / V0**2) / omega

    # Compute the RK4 steps
    k1_r = c2oomega * k * dt
    k1_k = 0.5 * omega * f(r[0], r[1], r[2], A, B, C, type) * dt

    k2_r = c2oomega * (k + 0.5 * k1_k) * dt
    k2_k = 0.5 * omega * f(r[0] + 0.5 * k1_r[0], r[1] + 0.5 * k1_r[1], r[2] + 0.5 * k1_r[2], A, B, C, type) * dt

    k3_r = c2oomega * (k + 0.5 * k2_k) * dt
    k3_k = 0.5 * omega * f(r[0] + 0.5 * k2_r[0], r[1] + 0.5 * k2_r[1], r[2] + 0.5 * k2_r[2], A, B, C, type) * dt

    k4_r = c2oomega * (k + k3_k) * dt
    k4_k = 0.5 * omega * f(r[0] + k3_r[0], r[1] + k3_r[1], r[2] + k3_r[2], A, B, C, type) * dt

    # Update the position and wave vector based on RK4 integration
    r_next = r + (k1_r + 2 * k2_r + 2 * k3_r + k4_r) / 6
    k_next = k + (k1_k + 2 * k2_k + 2 * k3_k + k4_k) / 6

    return r_next, k_next

def ray_trace_3d(r, k, dt, steps, A=1, B=0.01, C=0, type='exponential'):
    """
    Ray tracing in a 3D plasma density gradient with specified density distribution.
    r: Initial position of the ray (3D vector).
    k: Initial wave vector (direction) of the ray (3D vector).
    dt: Step size for integration.
    steps: Number of integration steps.
    A, B, C: Parameters for the plasma density distribution.
    type: The type of plasma density distribution, e.g., 'exponential'.
    """
    path = [r]
    
    for _ in range(steps):
        r, k = rk4_step_3d(compute_gradient_3d, r, k, dt, A, B, C, type)
        path.append(r)
        
    return np.array(path)

# Some physical parameters
elc = 1.602177e-19  # Elementary charge in Coulombs
m_e = 9.10938356e-31  # Electron mass in kg
mu_0 = 4 * np.pi * 1e-7  # Magnetic constant in Tm/A
clight = 3e8  # Speed of light in m/s
lambda_laser = 532e-9  # Wavelength in meters
omega_laser = 2 * np.pi * clight / lambda_laser  # Angular frequency (rad/s)

# Define the nondimensionalization 
L0 = lambda_laser  # Length scale in meters
L0 = 1e-3  # Length scale in meters
V0 = clight  # Time scale in seconds
T0 = L0/V0  # Velocity scale in m/s
n0 = find_critical_dens(lambda_laser*10**6)

# Define the plasma density parameters
I = 200e3 # current in Amperes
R = 85e-6 # radius in microns
# R = 150e-6 # radius in microns

Te = 1 # electron temperature in eV
Np = mu_0*I**2/(4*np.pi**2*R**2*elc*Te)  # Peak plasma density in 1/m^3

# Generate a 3D meshgrid for density plot
x = np.linspace(-1, 1, 500, dtype=np.float32)  
y = np.linspace(-1, 1, 500, dtype=np.float32)
z = np.linspace(-1, 1, 500, dtype=np.float32) 
X, Y, Z = np.meshgrid(x, y, z)  # This Z is different from plasma density Z, consider renaming for clarity
A = Np/n0
# A=1.0
B = R/L0
C = 0
type = 'bennett'

# Assuming plasma_density_3d is the 3D version of the plasma density function
# Density calculation in 3D would require adapting or iterating over 3D space
# Direct 3D visualization or analysis here would be complex and might require specific approaches

# Now none-dimensionalize the parameters
omega = omega_laser * T0  # Nondimensionalized angular frequency

# Compute the ray's path in 3D
dt = 0.0001  # Step size
steps = 10000  # Number of steps

# Initialize based on laser characteristics
fmec = 250e-3  # Focal length in meters
rodsize = 72e-3/2  # Radius of the rod in meters
angle = np.arctan(rodsize/fmec)  # Angle of the rod in radians
r50 = 68.33e-6/L0  # Radius of the minimum beam in meters
# r50 = 121.17e-6/L0  # Radius of the minimum beam in meters
dist = 1*B
ll = np.abs(x.min())-dist
dL = np.tan(angle)*ll

r0A = [np.array([r50+dL, -1, 0]), np.array([-r50-dL, -1, 0])]  # Starting point of the ray (x, y, z) already nondimensionalized

paths = []
plotter = None
for r0 in r0A:
    r0x = r0[0]
    if r0x > 0:
      direction = np.array([-dL, ll, 0])
    else:
      direction = np.array([dL, ll, 0])           

    # direction = np.array([B-r0[0], ll, 0])  # Initial direction, now in 3D
    k0 = L0 * direction / np.linalg.norm(direction) * omega_laser / clight  # Initial nondimensionalized direction, adjusted for 3D
    path = ray_trace_3d(r0, k0, dt, steps, A, B, C, type)
    paths.append(path)

    xs = path[:, 0]
    ys = path[:, 1]
    zs = path[:, 2]

    # Find minimum distance to the origin
    min_distance = np.min(np.linalg.norm(path, axis=1))
    print("Minimum distance to origin: ", min_distance)

    # Find location of the minimum distance
    min_distance_index = np.argmin(np.linalg.norm(path, axis=1))
    min_distance_location = path[min_distance_index]
    print("Location of minimum distance: ", min_distance_location)

    # Find the plasma density at this minimum location 
    min_distance_density = plasma_density_3d(path[min_distance_index, 0], path[min_distance_index, 1], path[min_distance_index, 2], A=A, B=B, C=C, type=type)

    # min_distance_density = plasma_density_3d(*min_distance_location, A, B, C, type)
    print("Plasma density at minimum distance: ", min_distance_density)

    # plot_ray_path_3d(path, plasma_density_3d=plasma_density_3d, A=A, B=B, C=C, type=type, X=x, Y=y, Z=z)
    plotter = plot_ray_path_3d_pyvis(plotter, path, plasma_density_3d=plasma_density_3d, A=A, B=B, C=C, type=type, x=x, y=y, z=z)
# Optionally, you can add axes grid with labels (this also adds labels)
# Add axis labels
# Show the axes grid with custom labels and other optional customizations
# plotter.show_grid()
# plotter.show_axes()
# plotter.view_isometric()
# plotter.set_position([2, -2, -4])
plotter.show_bounds(
    grid='back',
    location='outer',
    ticks='both',
    n_xlabels=2,
    n_ylabels=2,
    n_zlabels=2,
    xtitle='X',
    ytitle='Y',
    ztitle='Z',
)
plotter.reset_camera( bounds=(-1, 1, -1, 1, -1, 1))
plotter.camera_position = [(3, -3, 3), (0, 0, 0), (0, 0, 1)]

# plotter.view_vector([2, -2, 4])
plotter.remove_scalar_bar()
plotter.screenshot(r'D:\\'+'XSPL/Proposals/LCLS/2024/Figures\Plots/bennett_r_85um_with_mec_cpp_150um_nobar.png')  
plotter.save_graphic(r'D:\\'+'XSPL/Proposals/LCLS/2024/Figures\Plots/bennett_r_85um_with_mec_cpp_150um_nobar.svg') 
plotter.show()
