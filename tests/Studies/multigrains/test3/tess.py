import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

# Load the initial grain positions from file
# Each row is assumed to contain (x, y, z) coordinates
pts = np.loadtxt("initial-grain.txt")

# Project the 3D points onto the XY plane
# The z-coordinate is discarded
pts2d = pts[:, :2]

# Compute the 2D Voronoi tessellation
# using the projected (x, y) coordinates
vor = Voronoi(pts2d)

# Create a figure for visualization
fig, ax = plt.subplots(figsize=(8, 8))

# Plot the Voronoi diagram
# - show_vertices=False hides Voronoi vertices
# - line_width controls the edge thickness
voronoi_plot_2d(
    vor,
    ax=ax,
    show_vertices=False,
    line_width=1
)

# Plot the original seed points (Voronoi generators)
ax.scatter(
    pts2d[:, 0],
    pts2d[:, 1],
    c="red",
    s=15
)

# Use equal scaling on both axes
# so that distances are represented correctly
ax.set_aspect("equal")

# Display the figure
plt.show()