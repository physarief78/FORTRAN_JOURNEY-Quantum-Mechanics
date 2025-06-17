import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import os

# Simulation parameters (must match Fortran)
nt       = 10000
nsteps   = 10
nx, ny   = 1000, 1000
Lx, Ly   = 20e-9, 20e-9
dx, dy   = Lx/(nx-1), Ly/(ny-1)
dt       = 1.0e-18

# Physical constants for energy computation
ħ = 1.0545718e-34   # J·s
m = 9.10938356e-31  # kg (electron mass)
eV = 1.602176634e-19

# Gaussian-pillar lattice parameters (for overlay)
a_latt     = 0.5e-9      # 0.5 nm
x_lat_max  = 2.5e-9      # half the box in x
sigma_latt = 0.1e-9      # 0.1 nm

# Paths to the binary files
folder    = "output_data"
bin_psi   = os.path.join(folder, "wavepacket_psi_2d.bin")           # complex psi
bin_dens  = os.path.join(folder, "wavepacket_animation_2d.bin")    # density only


# --- Load density data (as before) ---
raw = np.fromfile(bin_dens, dtype=np.float64)
block_size = 1 + nx*ny
nframes = raw.size // block_size
if raw.size % block_size != 0:
    raise ValueError(f"Density file size {raw.size} is not a multiple of {block_size}.")

# Each frame = 1 time + nx*ny density values
block_size = 1 + nx*ny
nframes = raw.size // block_size
if raw.size % block_size != 0:
    raise ValueError(f"Binary file size {raw.size} is not a multiple of {block_size}.")

# Reshape into [nframes, 1+nx*ny]
data = raw.reshape(nframes, block_size)
times = data[:, 0]
dens_y = data[:, 1:].reshape(nframes, nx, ny)
density_data = (data[:, 1:].reshape(nframes, nx, ny)).transpose(0, 2, 1)

# Zero out the hard boundaries
density_data[:,  0, :] = 0
density_data[:, -1, :] = 0
density_data[:, :,  0] = 0
density_data[:, :, -1] = 0

# --- Load complex wavefunction psi_data ---
# Assuming psi file stores: for each frame, 2*nx*ny float64's = [Re ψ[0,0], Im ψ[0,0], Re ψ[1,0], Im ψ[1,0], ...]
# — Mem-map ψ file & prepare slicing along y=0 —
psi_block = 2*nx*ny
psi_mem   = np.memmap(bin_psi, dtype=np.float64, mode='r')
j_center  = np.argmin(np.abs(np.linspace(-Ly/2, Ly/2, ny)))  # index for y≈0

# Coordinate grids
x = np.linspace(-Lx/2, Lx/2, nx)
y = np.linspace(-Ly/2, Ly/2, ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Draw Gaussian‐pillar lattice
pillar_radius = sigma_latt

# Precompute lattice indices
nx_cells = int(x_lat_max / a_latt)
ny_cells = int((Ly/2)    / a_latt)

# --- 2D Snapshots at specified times ---
snap_times  = [times[0], 4.10e-16, 7.80e-16, 1.00e-15]
snap_labels = ["t_init",  "t4_10e16", "t7_80e16", "t1_00e15"]

for t_snap, label in zip(snap_times, snap_labels):
    k_snap = np.argmin(np.abs(times - t_snap))

    fig, ax = plt.subplots(figsize=(6,6))
    im = ax.imshow(
        density_data[k_snap].T,
        extent=(-Lx/2, Lx/2, -Ly/2, Ly/2),
        origin='lower',
        cmap='viridis',
        vmin=0, vmax=density_data.max()*0.1,
    )
    ax.set_title(f"Electron's Wavepacket Evolution @ t = {times[k_snap]:.2e} s")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

    # Draw lattice overlay
    for ix in range(nx_cells+1):
        cx = ix * a_latt
        if cx < -Lx/2 or cx > Lx/2:
            continue
        for iy in range(-ny_cells, ny_cells+1):
            cy = iy * a_latt
            if abs(cy) <= Ly/2:
                circ = Circle(
                    (cx, cy),
                    radius=sigma_latt,
                    facecolor='white',
                    edgecolor='black',
                    alpha=0.6,
                    linewidth=0.7
                )
                ax.add_patch(circ)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Probability Density")

    fname2d = f"wavepacket2d_{label}.png"
    plt.tight_layout()
    plt.savefig(fname2d, dpi=300)
    print(f"Saved 2D snapshot → {fname2d}")
    plt.close(fig)

# --- Full 2D Animation ---
fig_anim, ax_anim = plt.subplots(figsize=(6,6))
cax = ax_anim.imshow(
    density_data[0].T,
    extent=(-Lx/2, Lx/2, -Ly/2, Ly/2),
    origin='lower',
    cmap='viridis',
    vmin=0, vmax=density_data.max()*0.1,
)
cbar = fig_anim.colorbar(cax, ax=ax_anim)
cbar.set_label("Probability Density")

# redraw lattice once
for ix in range(nx_cells+1):
    cx = ix * a_latt
    if cx < -Lx/2 or cx > Lx/2:
        continue
    for iy in range(-ny_cells, ny_cells+1):
        cy = iy * a_latt
        if abs(cy) <= Ly/2:
            circ = Circle(
                (cx, cy),
                radius=sigma_latt,
                facecolor='white',
                edgecolor='black',
                alpha=0.6,
                linewidth=0.7
            )
            ax_anim.add_patch(circ)

ax_anim.set_xlim(-Lx/2, Lx/2)
ax_anim.set_ylim(-Ly/2, Ly/2)
ax_anim.set_xlabel("x (m)")
ax_anim.set_ylabel("y (m)")
ax_anim.set_title(f"Electron Scattering Simulation @ t = {times[0]:.2e} s")

def update(frame_idx):
    cax.set_data(density_data[frame_idx].T)
    ax_anim.set_title(f"Electron Scattering Simulation @ t = {times[frame_idx]:.2e} s")
    return (cax,)

ani = FuncAnimation(fig_anim, update, frames=nframes, blit=True, interval=20)
ani.save("wavepacket_2d.gif", fps=30)

plt.show()

psi_block = 2 * nx * ny
bytes_pf   = psi_block * 8

ħ = 1.0545718e-34
m = 9.10938356e-31
y_target = 1.0e-9
jy = np.argmin(np.abs(y - y_target))

targets = [
    (4.10e-16,  1.75e-9, "t4_10e16_x1_75nm"),
    (7.80e-16, -5.00e-9, "t7_80e16_x_minus5nm"),
    (1.00e-15,  7.50e-9, "t1_00e15_x7_5nm"),
]

colors = ['C0','C1','C2']

# Separate plots
# Separate plots with energy in eV
for (t_target, x_target, suffix), color in zip(targets, colors):
    # 1) find the frame and x‐index
    k   = np.argmin(np.abs(times - t_target))
    ix  = np.argmin(np.abs(x     - x_target))

    # 2) extract real density slice (for plotting)
    dens_1d = dens_y[k, :, ix]

    # 3) load the matching ψ‐frame
    with open(bin_psi, 'rb') as f:
        f.seek(k * bytes_pf)
        buf = np.fromfile(f, dtype=np.float64, count=psi_block)
    psi2d = (buf[0::2] + 1j*buf[1::2]).reshape((nx, ny), order='F')

    # 4) slice ψ along y at this x
    ψ1d = psi2d[ix, :]    # shape (ny,)

    # 5) compute second derivative in y
    d2ψ = np.empty_like(ψ1d)
    d2ψ[1:-1] = (ψ1d[2:] - 2*ψ1d[1:-1] + ψ1d[:-2]) / dy**2
    d2ψ[0]    = (ψ1d[1]   - 2*ψ1d[0]    + 0      ) / dy**2
    d2ψ[-1]   = (0        - 2*ψ1d[-1]  + ψ1d[-2]) / dy**2

    # 6) compute kinetic‐energy expectation ⟨E⟩ in eV
    num_J  = - (ħ**2/(2*m)) * np.trapezoid(np.conj(ψ1d) * d2ψ, y)
    den     = np.trapezoid(np.abs(ψ1d)**2, y)
    E_exp_J= (num_J/den).real
    E_exp_eV = E_exp_J / eV

    # 7) plot density vs y and annotate E in eV
    plt.figure(figsize=(6,4))
    plt.plot(y*1e9, dens_1d, color=color, lw=2)
    plt.xlabel("y (nm)")
    plt.ylabel("Probability Density")
    plt.title(
        "Diffraction Pattern of Electron's Wavepacket @\n"
        f"t = {times[k]:.2e} s, x = {x_target*1e9:.2f} nm\n"
        f"⟨E⟩ = {E_exp_eV:.2f} eV"
    )
    plt.grid(True)
    plt.tight_layout()

    fname = f"density_vs_y_{suffix}.png"
    plt.savefig(fname, dpi=300)
    print(f"Saved → {fname} (⟨E⟩ = {E_exp_eV:.2f} eV)")
    plt.close()

# Combined plot
plt.figure(figsize=(7,5))
for (t_target, x_target, suffix), color in zip(targets, colors):
    frame_idx = np.argmin(np.abs(times - t_target))
    ix        = np.argmin(np.abs(x     - x_target))
    dens_1d   = dens_y[frame_idx, :, ix]
    plt.plot(y*1e9, dens_1d, color=color, linewidth=2,
             label=f"t={times[frame_idx]:.2e}s, x={x_target*1e9:.2f}nm")

plt.xlabel("y (nm)")
plt.ylabel("Probability Density")
plt.title("Diffraction Pattern of Electron's Wavepacket @ Selected (t, x)")
plt.legend()
plt.grid(True)
plt.tight_layout()

combined_fname = "density_vs_y_combined.png"
plt.savefig(combined_fname, dpi=300)
print(f"Saved → {combined_fname}")
plt.show()

# y‐slice target (in meters) and find its grid index
y_target = 1.0e-9
jy_slice  = np.argmin(np.abs(y - y_target))

# region masks on x
maskA = x <  0
maskB = (x >= 0) & (x <= 2.5e-9)
maskC = x >  2.5e-9

# ensure output folder exists
os.makedirs(folder, exist_ok=True)

# Add initial time to targets
initial = (times[0], 0.0, "t_init")   # x_target is unused for y=1nm slice
targets_with_init = [initial] + targets
colors_with_init  = ['C3'] + colors

for (t_target, x_target, suffix), color in zip(targets_with_init, colors_with_init):
    # Find snapshot frame
    k = np.argmin(np.abs(times - t_target))

    # Read one ψ-frame
    with open(bin_psi, 'rb') as f:
        f.seek(k * bytes_pf)
        raw_frame = np.fromfile(f, dtype=np.float64, count=psi_block)

    # Reconstruct & transpose
    psi2d = (raw_frame[0::2] + 1j*raw_frame[1::2]).reshape((nx, ny), order='F').T

    # 1D slice at y ≃ 1 nm
    ψ1d   = psi2d[:, jy_slice]

    # Density and second derivative
    dens_x = np.abs(ψ1d)**2
    d2ψ    = np.empty_like(ψ1d)
    d2ψ[1:-1] = (ψ1d[2:] - 2*ψ1d[1:-1] + ψ1d[:-2]) / dx**2
    d2ψ[0]    = (ψ1d[1]   - 2*ψ1d[0]    + 0      ) / dx**2
    d2ψ[-1]   = (0        - 2*ψ1d[-1]  + ψ1d[-2]) / dx**2

    # Regional energies
    # compute regional energies in eV
    E = {}
    for name, mask in zip(['A','B','C'], [maskA,maskB,maskC]):
        num    = - (ħ**2/(2*m)) * np.trapezoid(np.conj(ψ1d[mask])*d2ψ[mask], x[mask])
        den    = np.trapezoid(dens_x[mask], x[mask])
        E[name] = (num/den).real / eV  # now in eV

    # plot
    plt.figure(figsize=(6,4))
    plt.plot(x*1e9, dens_x, color=color, lw=2)
    for xc in [0, 2.5e-9]:
        plt.axvline(x=xc*1e9, color='k', ls='--')

    ymax = dens_x.max()
    for name, mask, ypos in zip(['A','B','C'], [maskA,maskB,maskC], [0.9,0.8,0.7]):
        x_mid = x[mask].mean()*1e9
        plt.text(x_mid, ypos*ymax,
                 f"E{name} = {E[name]:.2f} eV",
                 color=color, ha='center', fontsize=10)

    plt.xlabel("x (nm)")
    plt.ylabel("Probability Density")
    plt.title(f"Electron's Wavepacket Evolution @ t={times[k]:.2e}s, y={y_target*1e9:.1f} nm")
    plt.grid(True)
    plt.tight_layout()

    outname = f"density_energy_{suffix}.png"
    plt.savefig(outname, dpi=300)
    print(f"Saved → {outname}")
    plt.close()

    
plt.figure(figsize=(7,5))
for (t_target, x_target, suffix), color in zip(targets, colors):
    frame_idx = np.argmin(np.abs(times - t_target))
    dens_x = dens_y[frame_idx, jy_slice, :]
    plt.plot(x*1e9, dens_x, color=color, lw=2,
             label=f"t={times[frame_idx]:.2e}s")
for xc in [0, 2.5e-9]:
        plt.axvline(x=xc*1e9, color='k', ls='--')
plt.xlabel("x (nm)")
plt.ylabel("Probability Density")
plt.title(f"Electron's Wavepacket Evolution Comparison\n@ different t with y = {y_target*1e9:.1f} nm")
plt.legend()
plt.grid(True)
plt.tight_layout()
combined_fname = "density_vs_x_y1nm_combined.png"
plt.savefig(combined_fname, dpi=300)
print(f"Saved → {combined_fname}")
plt.show()
