from proteus.iproteus import *
from proteus import Archiver
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle

import os
import errno


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


def post_process_case(case_name):
    archive = Archiver.XdmfArchive(".", case_name, readOnly=True)
    # Extract the number of time steps
    # Apparently Proteus creates a collection of grids
    M = 0
    for a in archive.tree.getroot().iter("Grid"):
        if a.attrib['GridType'] == "Collection":
            if a.attrib['Name'] == "Mesh_dgp2_Lagrange":
                Time = list(enumerate(a.iter("Time")))
                M = len(Time)

    # Obtain the position and velocities
    label_base = "/{0:s}{1:d}"

    x = read_from_hdf5(
        archive.hdfFile, label_base.format("nodesSpatial_Domain", 0))[0:-1, 0]
    u = read_from_hdf5(archive.hdfFile, label_base.format("u", M - 1))
    archive.close()
    return u[0:-1:3], x


def construct_mc_histogram(N_x=50, N_a=50):
    # Create the bins
    bin_bounds_x = np.linspace(0, 2 * np.pi, num=N_x)
    bin_bounds_a = np.linspace(-4, 4, num=N_a)
    X, A = np.meshgrid(bin_bounds_x, bin_bounds_a)
    H = np.zeros((N_a, N_x))

    # Get the solution
    N = 14000
    ns = range(N)
    for n in ns:
        file_name = "mc_run{0:05d}".format(n)
        print("Post-processing {0}...".format(file_name)),
        u_result, x_result = post_process_case(file_name)
        Hc, xedges, yedges = np.histogram2d(x_result, u_result, bins=[N_x, N_a], range=[[0, 2*np.pi], [-4,4]])
        H = H + Hc.T
        print 'done.'

    with open('mc_histogram.pickle', 'w') as f:
        pickle.dump([X, A, H], f)


construct_mc_histogram(100, 200)
