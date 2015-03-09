from proteus.iproteus import *
from proteus import Archiver
import pickle


def read_from_hdf5(hdfFile, label, dof_map=None):
    """
    Just grab the array stored in the node with label label and return it
    If dof_map is not none, use this to map values in the array
    If dof_map is not none, this determines shape of the output array
    """
    assert hdfFile is not None  # requires hdf5 for heavy data"
    vals = hdfFile.getNode(label).read()
    if dof_map is not None:
        dof = vals[dof_map]
    else:
        dof = vals

    return dof

archive = Archiver.XdmfArchive(".", "kinetic", readOnly=True)
label_base = "/{0:s}{1:d}"

M = 0
for a in archive.tree.getroot().iter("Grid"):
    if a.attrib['GridType'] == "Collection":
        if a.attrib['Name'] == "Mesh_dgp2_Lagrange":
            Time = list(enumerate(a.iter("Time")))
            M = len(Time)

p = read_from_hdf5(archive.hdfFile, label_base.format("p", M - 1))
x = read_from_hdf5(
    archive.hdfFile, label_base.format("nodes_dgp2_Lagrange", M - 1))

plot_histogram(x, p)
with open('xp.pickle', 'w') as f:
    pickle.dump((x, p), f)
