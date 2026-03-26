''' Verify that the temperatures at three given points match the expected values'''

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

# Read VTU file
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("Saves/Heat/Cycle000001/proc000000.vtu")
reader.Update()

# Get points and temperatures
data = reader.GetOutput()
points = vtk_to_numpy(data.GetPoints().GetData())
T_val = vtk_to_numpy(data.GetPointData().GetArray("T"))

# Define expected temperatures and points coordinates
# from J. P. Holman, Heat Transfer Tenth Edition. McGraw-Hill, pp. 111, Example 3-10 (2008). (corresponding to T4, T5 and T7 in example)
# and https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransferVerificationTests.html
points_ref = np.array([[0.005, 0.005, 0.0],[0.01,0.005,0.0],[0.005,0.0,0.0]])
T_ref = np.array([1092.37,1064.21,1111.38])

for i in range(3):
    id = np.where((points == points_ref[i]).all(axis=1))[0][0]

    # Check the temperature values against the expected values
    assert abs(T_val[id]-T_ref[i])/T_ref[i] < 0.01, f"Temperature at at least one point not matching expected value."