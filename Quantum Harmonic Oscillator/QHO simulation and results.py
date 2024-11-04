import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

# Define the filename and number of levels to plot
filename = 'energy_levels_dataV2'
num_levels = 5  # Ground state to fourth energy level

# Initialize lists to store data
x_values = None
probability_densities = []

# Read data from the file
with open(filename, 'r') as file:
    lines = file.readlines()
    energy_level_count = -1
    temp_prob_density = []

    for line in lines:
        if "Energy level" in line:
            # Append collected data for the previous level if exists
            if temp_prob_density:
                probability_densities.append(np.array(temp_prob_density))
                temp_prob_density = []
            energy_level_count += 1
            continue
        elif energy_level_count >= 0 and line.strip() and "Position" not in line:
            # Extract x and probability density values
            x, prob_density = map(float, line.split())
            if x_values is None:
                x_values = []
            x_values.append(x)
            temp_prob_density.append(prob_density)

    # Append the last set of data
    if temp_prob_density:
        probability_densities.append(np.array(temp_prob_density))

# Convert x_values to a numpy array and remove duplicates
x_values = np.unique(np.array(x_values))

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('black')   # Set the entire figure background to black
ax.set_facecolor('black')           # Set the plot (axes) background to black

# Get a colormap with distinct colors
cmap = get_cmap('rainbow', num_levels)

# Set custom labels for each energy level
energy_labels = ["Ground State", "First Excited State", "Second Excited State", "Third Excited State", "Fourth Excited State"]

for level in range(num_levels):
    ax.plot(x_values, probability_densities[level], label=energy_labels[level], color=cmap(level))

ax.set_xlim(-5, 2)                  # Set x-axis range from -5 to 2
ax.set_xlabel("Position (x)", color='white')     # Set x-axis label color to white
ax.set_ylabel("Probability Density", color='white')  # Set y-axis label color to white
ax.set_title("Probability Densities for Quantum Harmonic Oscillator Energy Levels", color='white')
ax.tick_params(colors='white')      # Set tick color to white

# Configure the legend with white text color for labels
legend = ax.legend(facecolor='black', edgecolor='white', labelcolor='white')
for text in legend.get_texts():
    text.set_color("white")         # Set legend text color to white

ax.grid(color='gray')               # Set grid color for better visibility on dark background

plt.show()








