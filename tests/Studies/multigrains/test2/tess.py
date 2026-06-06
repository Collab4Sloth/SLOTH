# from scipy.spatial import Voronoi
# import pyvista as pv
# import numpy as np

# pts = np.loadtxt("centroids.txt")

# vor = Voronoi(pts)

# plotter = pv.Plotter()

# plotter.add_points(
#     pts,
#     color='red',
#     point_size=10,
#     render_points_as_spheres=True
# )

# # Affiche les vertices Voronoï
# plotter.add_points(
#     vor.vertices,
#     color='blue',
#     point_size=5
# )

# plotter.show()


import numpy as np
import pyvista as pv
import pyvoro

# ==========================================================
# lecture des centroides
# ==========================================================

pts = np.loadtxt("centroids.txt")

# boite périodique 32x32x32
limits = [
    [0.0, 32.0],
    [0.0, 32.0],
    [0.0, 32.0]
]

# ==========================================================
# calcul Voronoi 3D
# ==========================================================

cells = pyvoro.compute_voronoi(
    pts,
    limits,
    4.0,
    periodic=[True, True, True]
)

# ==========================================================
# affichage
# ==========================================================

plotter = pv.Plotter()

# centroïdes rouges
plotter.add_points(
    pts,
    color="red",
    point_size=12,
    render_points_as_spheres=True
)

# cellules
for c in cells:

    vertices = np.array(c["vertices"])

    faces = []

    for face in c["faces"]:

        ids = face["vertices"]

        # format vtk :
        # [n,id1,id2,id3,...]
        faces.extend([len(ids)] + ids)

    poly = pv.PolyData(vertices, faces)

    plotter.add_mesh(
        poly,
        opacity=0.25,
        show_edges=True
    )

plotter.show()


# import numpy as np
# import pyvista as pv

# pts = np.loadtxt("centroids.txt")

# cloud = pv.PolyData(pts)

# plotter = pv.Plotter()
# plotter.add_mesh(cloud,
#                  render_points_as_spheres=True,
#                  point_size=12,
#                  color="red")

# plotter.show_grid()
# plotter.show()