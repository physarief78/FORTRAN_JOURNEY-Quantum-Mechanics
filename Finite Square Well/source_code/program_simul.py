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

# Define the finite square well potential
V0 = 100.0  # Height of the finite square well (same as in the Fortran code)
a = 1.0    # Half-width of the finite square well (same as in the Fortran code)
V = np.zeros_like(x)
V[np.abs(x) > a] = V0  # Potential is V0 outside the well and 0 inside

# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 8))
line_density, = ax.plot(x, densities[0, :], color="blue", lw=2, label="Probability Density")
line_real, = ax.plot(x, real_parts[0, :], color="green", lw=2, label="Real Part")
line_imag, = ax.plot(x, imag_parts[0, :], color="purple", lw=2, label="Imaginary Part")
line_potential, = ax.plot(x, V / np.max(V) * np.max(densities), color="red", lw=2, label="Finite Square Well Potential")
ax.set_xlim(-L/2, L/2)
ax.set_ylim(-2, 1.1 * np.max(densities))
ax.set_title("Wave Packet Dynamics in a Finite Square Well")
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
ani.save("wavepacket_finite_square_well.gif", fps=30, dpi=150)

# Show the animation
plt.show()