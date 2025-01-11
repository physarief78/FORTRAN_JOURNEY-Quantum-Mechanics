import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the data
filename = "wavepacket_harmonic.dat"
data = np.loadtxt(filename)

# Extract time steps and density values
time_steps = data[:, 0]
densities = data[:, 1:]

# Number of spatial grid points (nx) and time steps (nt)
nx = densities.shape[1]
nt = len(time_steps)

# Define the spatial grid
L = 3.0  # Length of the box (same as in the Fortran code)
x = np.linspace(-L/2, L/2, nx)

# Define the harmonic potential
m = 1.0  # Mass (same as in the Fortran code)
omega = 10.0  # Angular frequency of the harmonic oscillator
V = 0.5 * m * omega**2 * x**2  # Harmonic potential

# Create a figure and axis
fig, ax = plt.subplots()
line_density, = ax.plot(x, densities[0, :], color="blue", lw=2, label="Probability Density")
line_potential, = ax.plot(x, V / np.max(V) * np.max(densities), color="red", lw=2, label="Harmonic Potential")
ax.set_xlim(-L + 1, L - 1)
ax.set_ylim(-1, 1.1 * np.max(densities))
ax.set_title("Wave Packet Dynamics in a Harmonic Potential")
ax.set_xlabel("Position (x)")
ax.set_ylabel("Probability Density")
ax.legend()

# Update function for the animation
def update(frame):
    line_density.set_ydata(densities[frame, :])
    ax.set_title(f"Wave Packet Dynamics at t = {time_steps[frame]:.3f}")
    return line_density, line_potential

# Create the animation
ani = FuncAnimation(fig, update, frames=nt, interval=30, blit=True)

# Save the animation as a video or GIF (optional)
ani.save("wavepacket_harmonic_with_potential.gif", fps=30, dpi=150)

# Show the animation
plt.show()
