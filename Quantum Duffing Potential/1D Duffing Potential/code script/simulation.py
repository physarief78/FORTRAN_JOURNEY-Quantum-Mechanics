import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the data
filename = "wavepacket_duffing.dat"  # Updated file for Duffing potential output
data = np.loadtxt(filename)

# Extract time steps and density values
time_steps = data[:, 0]
densities = data[:, 1:]

# Number of spatial grid points (nx) and time steps (nt)
nx = densities.shape[1]
nt = len(time_steps)

# Define the spatial grid
L = 3.0  # Length of the box (same as in the Fortran code)
x = np.linspace(-L, L, nx)

# Define the Duffing potential
alpha = -10.0  # Duffing potential quadratic term coefficient (same as in Fortran)
beta = 2.0    # Duffing potential quartic term coefficient (same as in Fortran)
V = 2.0 +  (alpha * x**2 + beta * x**4)  # Duffing potential

# Normalize the potential for better visualization alongside the density
V_normalized = V / np.max(V) * np.max(densities)

# Create a figure and axis
fig, ax = plt.subplots()
line_density, = ax.plot(x, densities[0, :], color="blue", lw=2, label="Probability Density")
line_potential, = ax.plot(x, V_normalized, color="red", lw=2, label="Duffing Potential")
ax.set_xlim(-L, L)
ax.set_ylim(-1, 1.1 * np.max(densities))
ax.set_title("Wave Packet Dynamics in a Duffing Potential")
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
ani.save("wavepacket_duffing_with_potential.gif", fps=30, dpi=150)

# Show the animation
plt.show()
