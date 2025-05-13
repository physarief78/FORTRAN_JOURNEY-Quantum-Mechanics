import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

# Load the data from the Fortran output file
filename = "output_data/wavepacket_animation_2d.dat"
raw = np.loadtxt(filename)

# Simulation parameters (must match Fortran)
nt       = 2000
nx, ny   = 200, 200
Lx, Ly   = 2.5, 2.5
dx, dy   = Lx/nx, Ly/ny
dt       = 0.00001
nsteps   = 2

# Determine how many frames were saved
nframes = raw.shape[0]  # should be nt/nsteps = 500

# Extract times and reshape density
times        = raw[:, 0]
density_data = raw[:, 1:].reshape(nframes, nx, ny)

# Zero out hard boundaries
density_data[:,  0, :] = 0
density_data[:, -1, :] = 0
density_data[:, :,  0] = 0
density_data[:, :, -1] = 0

# Gaussian-pillar lattice parameters (must match Fortran)
V0         = 5000.0     # pillar height (not used in plotting)
a_latt     = 0.15       # lattice constant
x_lat_max  = 0.75       # lattice extends from x=0 to x=0.75
sigma_latt = 0.03       # width of each Gaussian pillar

# Compute how many cells in x and y
nx_cells = int(x_lat_max / a_latt)
ny_cells = int((Ly/2)     / a_latt)

# Create coordinate grids for density plotting
x = np.linspace(-Lx/2, Lx/2, nx)
y = np.linspace(-Ly/2, Ly/2, ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Set up the figure
fig, ax = plt.subplots(figsize=(6,6))

# Plot the first density frame
cax = ax.imshow(
    density_data[0].T,
    extent=(-Lx/2, Lx/2, -Ly/2, Ly/2),
    origin='lower',
    cmap='viridis',
    vmin=0, vmax=density_data.max()*0.1,
    alpha=1.0
)
cbar = fig.colorbar(cax, ax=ax)
cbar.set_label("Probability Density")

# Overlay each pillar as a translucent circle
pillar_radius = 1 * sigma_latt
for ix in range(nx_cells+1):
    cx = ix * a_latt
    if cx < 0 or cx > x_lat_max: 
        continue
    for iy in range(-ny_cells, ny_cells+1):
        cy = iy * a_latt
        # Only draw circles that lie inside the simulation box
        if abs(cy) <= Ly/2:
            circ = Circle(
                (cx, cy), 
                radius=pillar_radius, 
                color='white', 
                alpha=0.6, 
                edgecolor='black', 
                linewidth=0.7
            )
            ax.add_patch(circ)

ax.set_xlim(-Lx/2, Lx/2)
ax.set_ylim(-Ly/2, Ly/2)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("t = 0.00000 s")

# Animation update function
def update(i):
    cax.set_data(density_data[i].T)
    ax.set_title(f"t = {times[i]:.5f} s")
    return (cax,)

# Create and save the animation
ani = FuncAnimation(fig, update, frames=nframes, blit=True, interval=20)
ani.save("wavepacket_2d_bragg_diffraction_bot.gif", fps=30)

plt.show()
