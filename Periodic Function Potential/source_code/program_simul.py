import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the data
filename_density = "wavepacket_animation.dat"
filename_real = "real_parts_wavefunction.dat"
filename_imag = "imag_parts_wavefunction.dat"

data_density = np.loadtxt(filename_density)
data_real = np.loadtxt(filename_real)
data_imag = np.loadtxt(filename_imag)

# Extract time steps and values
time_steps = data_density[:, 0]
densities = data_density[:, 1:]
real_parts = data_real[:, 1:]
imag_parts = data_imag[:, 1:]

# Number of spatial grid points (nx) and time steps (nt)
nx = densities.shape[1]
nt = len(time_steps)

# Define the spatial grid
L = 10.0  # Length of the box (same as in the Fortran code)
x = np.linspace(-L/2, L/2, nx)

# Define the periodic potential
V0 = 0.25  # Amplitude of the periodic potential (same as in the Fortran code)
a = 1.0   # Periodicity of the potential (same as in the Fortran code)
V = V0 * np.cos(2.0 * np.pi * x / a) + V0 # Periodic potential

# Adjust the figure size here
fig_width = 12  # Width of the figure in inches
fig_height = 8  # Height of the figure in inches
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plot the initial wavefunction and potential
line_density, = ax.plot(x, densities[0, :], color="blue", lw=2, label="Probability Density")
line_real, = ax.plot(x, real_parts[0, :], color="green", lw=2, label="Real Part")
line_imag, = ax.plot(x, imag_parts[0, :], color="purple", lw=2, label="Imaginary Part")
line_potential, = ax.plot(x, V, color="red", lw=2, label="Periodic Potential")

# Set plot limits and labels
ax.set_xlim(-L/2, L/2)
ax.set_ylim(-2, 1.1 * np.max(densities))
ax.set_title("Wave Packet Dynamics in a Periodic Potential")
ax.set_xlabel("Position (x)")
ax.set_ylabel("Amplitude")
ax.legend()

# Update function for the animation
def update(frame):
    line_density.set_ydata(densities[frame, :])
    line_real.set_ydata(real_parts[frame, :])
    line_imag.set_ydata(imag_parts[frame, :])
    ax.set_title(f"Wave Packet Dynamics at t = {time_steps[frame]:.3f}")
    return line_density, line_real, line_imag, line_potential

# Create the animation
ani = FuncAnimation(fig, update, frames=nt, interval=30, blit=True)

# Save the animation as a video or GIF (optional)
ani.save("wavepacket_periodic_potential.gif", fps=30, dpi=150)

# Show the animation
plt.show()