import gmsh
import math
import sys
gmsh.initialize()
gmsh.model.add("2D_mesh")
geom = gmsh.model.occ


radius=0.00465
angle = math.pi/4.
h=8.e-5
hh=2.e-5
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 1)  
gmsh.model.geo.addPoint(radius, 0.0, 0.0, h, 2) 
gmsh.model.geo.addPoint(radius*math.cos(angle), radius*math.sin(angle),0.0, h, 3) 


gmsh.model.geo.addLine(1, 2, 4)
gmsh.model.geo.addLine(1, 3, 5)
gmsh.model.geo.addCircleArc(2, 1, 3, 6)

gmsh.model.geo.addCurveLoop([5,-6,-4],7) 
# gmsh.model.geo.addPlaneSurface([7], 8)

# Inclusion
gmsh.model.geo.addPoint(0.002, 0.0008, 0.0, hh, 9)  # center
gmsh.model.geo.addPoint(0.0025, 0.0008, 0.0, hh, 11)  # right
gmsh.model.geo.addPoint(0.002, 0.0013, 0.0, hh, 12)  # top
gmsh.model.geo.addPoint(0.002, 0.0003, 0.0, hh, 13)  # bottom
gmsh.model.geo.addPoint(0.0015, 0.0008, 0.0, hh, 14)  # left
gmsh.model.geo.addCircleArc(11, 9, 13, 15)
gmsh.model.geo.addCircleArc(13, 9, 14, 16)
gmsh.model.geo.addCircleArc(14, 9, 12, 17)
gmsh.model.geo.addCircleArc(12, 9, 11, 18)
gmsh.model.geo.addCurveLoop([15,16,17,18], 19)


# Surfaces
gmsh.model.geo.addPlaneSurface([7, 19], 20)
gmsh.model.geo.addPlaneSurface([19], 21)


gmsh.model.geo.synchronize()

# dimension, [entities], tag
gmsh.model.addPhysicalGroup(2, [20],22)
gmsh.model.setPhysicalName(2, 22, "pellet")
gmsh.model.addPhysicalGroup(2, [21],23)
gmsh.model.setPhysicalName(2, 23, "cluster")
# 
gmsh.model.addPhysicalGroup(1, [4],1)
gmsh.model.setPhysicalName(1, 1, "bottom")
gmsh.model.addPhysicalGroup(1, [5],2)
gmsh.model.setPhysicalName(1, 2, "top")
gmsh.model.addPhysicalGroup(1, [6],3)
gmsh.model.setPhysicalName(1, 3, "external")

gmsh.model.mesh.generate(2)

# Preview mesh.
gmsh.fltk.run()

# Clear mesh and close gmsh API.
gmsh.clear()
gmsh.finalize()