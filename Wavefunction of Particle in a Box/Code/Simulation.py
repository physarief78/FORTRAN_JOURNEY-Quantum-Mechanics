import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the data from the Fortran output file
data = np.loadtxt("wavepacket_animation.dat")
time_steps = data.shape[0]  # Number of time steps
nx = data.shape[1] - 1      # Number of spatial points

# Extract time and density data
time = data[:, 0]
density = data[:, 1:].reshape((time_steps, nx))

# Set up the spatial grid (assuming symmetric box centered at 0)
L = 10.0  # Length of the box
x = np.linspace(-L/2, L/2, nx)

# Create the figure and axis
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(-L/2, L/2)
ax.set_ylim(0, np.max(density))
ax.set_xlabel("Position x")
ax.set_ylabel("Probability Density")
ax.set_title("Gaussian Wavepacket Evolution")

def init():
    """Initialize the animation frame."""
    line.set_data([], [])
    return line,

def update(frame):
    """Update the animation for each frame."""
    line.set_data(x, density[frame])
    ax.set_title(f"Gaussian Wavepacket Evolution (Time: {time[frame]:.2f})")
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=time_steps, init_func=init, blit=True)

# Save or display the animation
ani.save("wavepacket_evolution.gif", fps=30, writer="pillow")
# Alternatively, show the animation interactively:
plt.show()
