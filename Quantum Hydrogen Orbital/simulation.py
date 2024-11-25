import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

# Function to load data from file
def load_orbital_data(filename):
    data = np.loadtxt(filename)
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    z = data[:, 2].reshape(len(x), len(y))
    return x, y, z

# Load all data files for n=4, l=0 to 3, m=0 to l
data_files = sorted(glob.glob("orbital_data_4_*_*_normalized.dat"))

# Create output directory for images
output_dir = "orbital_images"
os.makedirs(output_dir, exist_ok=True)

# Initialize plot
fig, ax = plt.subplots(figsize=(8, 6))
x, y, z = load_orbital_data(data_files[0])  # Use the first file to initialize
X, Y = np.meshgrid(x, y)
contour = ax.contourf(X, Y, z, levels=50, cmap="seismic_r") 
colorbar = plt.colorbar(contour, ax=ax, label="Wavefunction Amplitude (ψ)")

# Add a title to update separately
title = ax.set_title("")

# Update function for animation and saving
def save_frame(frame):
    global contour
    for c in contour.collections:
        c.remove()  # Clear previous contours

    # Load data for the current frame
    x, y, z = load_orbital_data(data_files[frame])
    contour = ax.contourf(X, Y, z, levels=50, cmap="seismic_r") 
    
    # Extract quantum numbers from filename
    filename = data_files[frame]
    # Filename format: orbital_data_4_<l>_<m>_normalized.dat
    parts = filename.split("_")  # Split filename on '_'
    l = int(parts[3])  # Extract l (corrected index)
    m = int(parts[4])  # Extract m (corrected index)

    # Update title dynamically
    title.set_text(f"Hydrogen Quantum State of n=4, l={l}, m={m}")
    
    # Save the current frame as an image
    output_filename = os.path.join(output_dir, f"orbital_n4_l{l}_m{m}.png")
    plt.savefig(output_filename, dpi=300)
    print(f"Saved: {output_filename}")

# Save all frames as individual images
for frame in range(len(data_files)):
    save_frame(frame)

print(f"All frames saved in '{output_dir}' directory.")

# Update function for animation
def update(frame):
    global contour
    for c in contour.collections:
        c.remove()  # Clear previous contours

    # Load data for the current frame
    x, y, z = load_orbital_data(data_files[frame])
    contour = ax.contourf(X, Y, z, levels=50, cmap="seismic_r") 
    
    # Extract quantum numbers from filename
    filename = data_files[frame]
    # Filename format: orbital_data_4_<l>_<m>_normalized.dat
    parts = filename.split("_")  # Split filename on '_'
    l = int(parts[3])  # Extract l (corrected index)
    m = int(parts[4])  # Extract m (corrected index)

    # Update title dynamically
    title.set_text(f"Hydrogen Quantum State of n=4, l={l}, m={m}")
    return contour.collections + [title]

# Create animation
ani = FuncAnimation(fig, update, frames=len(data_files), blit=False, interval=1000)

# Display or save animation
ani.save("hydrogen_orbital_animation.gif")  # Uncomment to save animation

plt.show()
