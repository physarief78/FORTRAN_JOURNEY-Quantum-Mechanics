import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Load data
data = np.loadtxt('wavepacket_animation_qho.dat')
time_steps = data[:, 0]
nx, ny = 100, 100  # Grid size (same as in Fortran code)
probability_density = data[:, 1:].reshape(len(time_steps), nx, ny)

# Define the spatial grid
Lx, Ly = 3.0, 3.0  # Same as in Fortran code
x = np.linspace(-Lx / 2, Lx / 2, nx)
y = np.linspace(-Ly / 2, Ly / 2, ny)
X, Y = np.meshgrid(x, y)

# Define the harmonic potential
m, omega = 1.0, 30.0  # Mass and angular frequency (same as in Fortran code)
harmonic_potential = 0.5 * m * omega**2 * (X**2 + Y**2)

# Normalize the probability density to match the scale of the harmonic potential
probability_density = probability_density / np.max(probability_density) * np.max(harmonic_potential)

# Create the figure and 3D axis
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Set the initial viewing angle
ax.view_init(elev=20, azim=20)  # Example: 30° elevation and 45° azimuth

# Initialize the plot with the potential and the first frame of the density
surf_density = ax.plot_surface(X, Y, probability_density[0], cmap='viridis', edgecolor='none', alpha=1.0, label='Probability Density')
surf_potential = ax.plot_surface(X, Y, harmonic_potential, cmap='inferno', edgecolor='none', alpha=0.3, label='Harmonic Potential')

# Set axis limits and labels
ax.set_title('Time Evolution of Quantum Harmonic Oscillator')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Z-axis')
ax.set_zlim(0, max(np.max(probability_density), np.max(harmonic_potential)) * 1.1)

# Add colorbars for probability density
fig.colorbar(surf_density, ax=ax, shrink=0.5, aspect=10, label='Probability Density')

# Update function for animation
def update(frame):
    global surf_density
    # Remove the previous probability density surface plot
    surf_density.remove()
    # Add a new surface plot for the current frame
    surf_density = ax.plot_surface(X, Y, probability_density[frame], cmap='viridis', edgecolor='none', alpha=1.0)
    ax.set_title(f'Time Evolution at t = {time_steps[frame]:.3f}')
    return [surf_density]

# Create the animation
ani = FuncAnimation(fig, update, frames=len(time_steps), interval=50, blit=False)

# Save the animation as a gif (optional)
ani.save('quantum_harmonic_oscillator_3d.gif', writer='pillow', fps=30)

# Show the animation
plt.show()


