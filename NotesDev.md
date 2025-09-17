# Development notes for integrating AMR in SLOTH

First step : ensure non conforming mesh -> modification of Spatial/Spatial.hpp  file
from ex15 in mfem
  EnsureNCMesh(true);
  // Make sure tet-only meshes are marked for local refinement.
  Finalize(true);



Work on AllenCahn/2D/test1  case 
Question about mesh coming from GMSH : do we specify in GMSH the NC ?