import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the data from the Fortran output file
filename = "wavepacket_animation_2d.dat"
data = np.loadtxt(filename)

# Parameters
nt = 1500  # Number of time steps
nx = 100  # Number of spatial grid points in x
ny = 100  # Number of spatial grid points in y
Lx, Ly = 2.0, 2.0  # Domain size
dx, dy = Lx / nx, Ly / ny

# Reshape the data into a 3D array (time, x, y)
density_data = data[:, 1:].reshape(nt, nx, ny)

# Create a mask for the infinite potential boundary (set values to zero at edges)
density_data[:, 0, :] = 0  # Top edge
density_data[:, -1, :] = 0  # Bottom edge
density_data[:, :, 0] = 0  # Left edge
density_data[:, :, -1] = 0  # Right edge

# Parameters for the double-slit potential
slit_width = 0.1  # Width of each slit
slit_separation = 0.2  # Separation between the centers of the slits

# Define the positions of the walls and slits
y = np.linspace(-Ly / 2, Ly / 2, ny)
x_wall_left = -0.01  # Slightly to the left of the slit center
x_wall_right = 0.01  # Slightly to the right of the slit center

# Define the slit regions
slit1_top = slit_separation / 4 + slit_width / 8
slit1_bottom = slit_separation / 4 - slit_width / 8
slit2_top = -slit_separation / 4 + slit_width / 8
slit2_bottom = -slit_separation / 4 - slit_width / 8

# Create a figure and axis for the animation
fig, ax = plt.subplots(figsize=(6, 6))
cax = ax.imshow(
    density_data[0], extent=(-Lx / 2, Lx / 2, -Ly / 2, Ly / 2), 
    cmap="viridis", origin="lower"
)
cbar = fig.colorbar(cax, ax=ax)
cbar.set_label("Probability Density")

# Overlay the walls using fill_betweenx
ax.fill_betweenx(
    y, x_wall_left, x_wall_right, 
    where=((y > slit1_top) | (y < slit1_bottom) & (y > slit2_top) | (y < slit2_bottom)), 
    color="white", alpha=1.0
)

ax.set_title("2D Gaussian Wavepacket Evolution with Double Slits")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Update function for animation
def update(frame):
    cax.set_data(density_data[frame])
    ax.set_title(f"2D Gaussian Wavepacket Evolution (t = {frame * 0.001:.4f}s)")
    return cax,

# Create the animation
ani = FuncAnimation(fig, update, frames=nt, blit=True)

# Save the animation as a file (optional)
ani.save("wavepacket_2d_animation.gif", fps=30)

# Show the animation
plt.show()



