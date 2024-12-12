import matplotlib.pyplot as plt
import numpy as np
import sys
from collections import Counter

# Function to parse input from stdin
def read_input():
    spike_angles = []
    intensities = []

    for line in sys.stdin:
        # Each line has the format "<angle> <intensity>"
        parts = line.strip().split()
        if len(parts) == 2:
            angle = float(parts[0])
            intensity = float(parts[1])
            spike_angles.append(angle)
            intensities.append(intensity)

    return spike_angles, intensities

# Function to plot a generic polygonal aperture (with a given number of sides)
def plot_polygon(center, size, num_sides, edge_intensities):
    # Create the angles for each vertex of the polygon
    theta = np.linspace(0, 2 * np.pi, num_sides + 1)  # num_sides+1 to close the polygon
    x = center[0] + size * np.cos(theta)
    y = center[1] + size * np.sin(theta)

    # Plot edges with intensity weighting
    for i in range(num_sides):
        intensity = edge_intensities[i % len(edge_intensities)] / max(edge_intensities)
        plt.plot(x[i:i + 2], y[i:i + 2], color=(1 - intensity, intensity, 0), linewidth=2)

# Function to determine the number of sides based on the spike angles
def determine_num_sides(spike_angles):
    # Sort spike angles
    spike_angles = sorted(spike_angles)

    # Calculate angular differences between consecutive angles
    angular_diffs = [(spike_angles[i] - spike_angles[i - 1]) % 360 for i in range(1, len(spike_angles))]

    # Find the most common angular difference (the most frequent spike spacing)
    common_spacing = Counter(angular_diffs).most_common(1)[0][0]

    # Calculate the number of sides
    num_sides = round(360 / common_spacing)
    return num_sides

# Main function to process data and plot
def main():
    # Read spike angles and intensities from C++ input via stdin
    spike_angles, intensities = read_input()

    if len(spike_angles) == 0 or len(intensities) == 0:
        print("No data received from C++ program.")
        return

    # Calculate edge directions (perpendicular to the spike direction)
    edge_angles = [(angle + 90) % 360 for angle in spike_angles]

    # Normalize intensities for coloring
    normalized_intensities = [intensity / max(intensities) for intensity in intensities]

    # Automatically determine the number of sides of the aperture based on spike angles
    num_sides = determine_num_sides(spike_angles)

    # Plot the reconstructed aperture
    plt.figure(figsize=(6, 6))
    plot_polygon(center=(0, 0), size=1, num_sides=num_sides, edge_intensities=normalized_intensities)
    plt.title(f"Reconstructed {num_sides}-Sided Aperture")
    plt.axis('equal')
    plt.show()

# Run the main function
if __name__ == "__main__":
    main()
