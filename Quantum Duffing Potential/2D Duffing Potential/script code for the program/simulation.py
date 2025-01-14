import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Load data
data = np.loadtxt('wavepacket_animation_duffing.dat')  # Updated file name
time_steps = data[:, 0]
nx, ny = 100, 100  # Grid size (same as in Fortran code)
probability_density = data[:, 1:].reshape(len(time_steps), nx, ny)

# Define the spatial grid
Lx, Ly = 3.0, 3.0  # Same as in Fortran code
x = np.linspace(-Lx / 2, Lx / 2, nx)
y = np.linspace(-Ly / 2, Ly / 2, ny)
X, Y = np.meshgrid(x, y)

# Define the Duffing potential
alpha, beta = -30.0, 25.0  # Same as in Fortran code
m, omega = 1.0, 30.0
duffing_potential = 10.0 + (alpha * (X**2 + Y**2) + beta * (X**4 + Y**4))

# Normalize the probability density to match the scale of the Duffing potential
probability_density = probability_density / np.max(probability_density) * np.max(duffing_potential)

# Create the figure and 3D axis
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Set the initial viewing angle
ax.view_init(elev=10, azim=45)  # Example: 20° elevation and 20° azimuth

# Initialize the plot with the potential and the first frame of the density
surf_density = ax.plot_surface(X, Y, probability_density[0], cmap='viridis', edgecolor='none', alpha=1.0, label='Probability Density')
surf_potential = ax.plot_surface(X, Y, duffing_potential, cmap='spring', edgecolor='none', alpha=0.3, label='Duffing Potential')

# Set axis limits and labels
ax.set_title('Time Evolution in Duffing Potential')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Z-axis')
ax.set_zlim(0, max(np.max(probability_density), np.max(duffing_potential)) * 1.1)

# Add colorbars for probability density
fig.colorbar(surf_density, ax=ax, shrink=0.5, aspect=10, label='Probability Density')

# Update function for animation
def update(frame):
    global surf_density
    # Remove the previous probability density surface plot
    surf_density.remove()
    # Add a new surface plot for the current frame
    surf_density = ax.plot_surface(X, Y, probability_density[frame], cmap='viridis', edgecolor='none', alpha=1.0)
    ax.view_init(elev=10, azim=2/3 * frame)
    ax.set_title(f'Time Evolution at t = {time_steps[frame]:.3f}')
    return [surf_density]

# Create the animation
ani = FuncAnimation(fig, update, frames=len(time_steps), interval=50, blit=False)

# Save the animation as a gif (optional)
ani.save('duffing_potential_wavepacket_3d_rotate.gif', writer='pillow', fps=30)

# Show the animation
plt.show()
