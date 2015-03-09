from proteus.iproteus import *
from proteus import Archiver
import numpy as np
from aux import set_plt_params
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib.patches import ConnectionPatch

import os
import errno


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


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


def plot_snapshot_on_axes(archive, axes, n=0, t=0):
    label_base = "/{0:s}{1:d}"

    x = read_from_hdf5(
        archive.hdfFile, label_base.format("nodes_dgp2_Lagrange", 0))[:, 0]
    u = read_from_hdf5(archive.hdfFile, label_base.format("u", n))
    axes.set_xlim([0, 2*np.pi])
    axes.set_ylim([-1.5, 4])
    axes.set_frame_on(False)
    axes.axes.get_yaxis().set_visible(False)
    axes.axes.get_xaxis().set_visible(False)
    axes.add_artist(mpl.lines.Line2D(
        (0, 2*np.pi), (0, 0), color='black', linewidth=0.5, alpha=0.5))
    axes.text(0.1, 3.0, "t={0:.1f}".format(float(t)))
    # Why third column is so smooth? What is the meaning of the columns?
    axes.plot(x, u, "k")


def make_solution_map_with_snapshots(archive, T=2, filename="map.pdf"):
    # Extract the number of time steps
    # Apparently Proteus creates a collection of grids
    M = 0
    for a in archive.tree.getroot().iter("Grid"):
        if a.attrib['GridType'] == "Collection":
            if a.attrib['Name'] == "Mesh_dgp2_Lagrange":
                Time = list(enumerate(a.iter("Time")))
                M = len(Time)

    # Number of nodes
    N = 1000
    temp = np.zeros((M * N, 1))

    label_base = "/{0:s}{1:d}"
    for ind in range(M):
        temp[ind*N:(ind+1)*N, 0] = read_from_hdf5(
            archive.hdfFile, label_base.format("u", ind))[0:-1:3]

    grid = temp.reshape((M, N))
    grid = grid[::-1, :]

    set_plt_params(relative_fig_width=1)
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, height_ratios=(1, 10))
    ax1 = plt.subplot(gs[1, 0])
    ax2 = plt.subplot(gs[0, 0])
    ax3 = plt.subplot(gs[:, 1])

    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$t$")
    ax1.set_xticks((0, 0.5, 1))
    ax1.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax1.set_yticks((0, 0.5, 1))
    ax1.set_yticklabels(("0", T/2, T))
    im = ax1.imshow(
        grid, extent=(0, 1, 0, 1), interpolation='nearest', cmap=cm.coolwarm,
        vmin=-1, vmax=2)

    cbar = plt.colorbar(im, orientation='horizontal', cax=ax2)
    cbar.solids.set_edgecolor("face")

    # I won't this area to be blank as I will put custom subplot here
    ax3.axis('off')
    ax4 = fig.add_axes([0.55, 0.71, 0.4, 0.2])
    ax5 = fig.add_axes([0.55, 0.38, 0.4, 0.2])
    ax6 = fig.add_axes([0.55, 0.05, 0.4, 0.2])
    t = [0, 1, 2]
    axes = [ax6, ax5, ax4]
    ind = 0
    for a in archive.tree.getroot().iter("Grid"):
        if a.attrib['GridType'] == "Collection":
            if a.attrib['Name'] == "Mesh_dgp2_Lagrange":
                for n, Time in enumerate(a.iter("Time")):
                    if ind >= len(t):
                        break
                    elif t[ind] <= float(Time.attrib['Value']):
                        print "Making snapshot for t = ", Time.attrib['Value']
                        plot_snapshot_on_axes(
                            archive, axes[ind], n=n, t=Time.attrib['Value'])
                        ind = ind + 1

    # Finishing touch - add arrows
    con = ConnectionPatch(
        xyA=(1, 0), xyB=(0, 0), coordsA="data", coordsB="data",
        axesA=ax1, axesB=ax6, arrowstyle="->", shrinkB=5, shrinkA=5)
    ax1.add_artist(con)
    con = ConnectionPatch(
        xyA=(1, 0.5), xyB=(0, 0), coordsA="data", coordsB="data",
        axesA=ax1, axesB=ax5, arrowstyle="->", shrinkB=5, shrinkA=5)
    ax1.add_artist(con)
    con = ConnectionPatch(
        xyA=(1, 1.0), xyB=(0, 0), coordsA="data", coordsB="data",
        axesA=ax1, axesB=ax4, arrowstyle="->", shrinkB=5, shrinkA=5)
    ax1.add_artist(con)
    fig.savefig(filename, bbox_inches='tight', transparent=True)


def post_process_case(case_name):
    archive = Archiver.XdmfArchive(".", case_name, readOnly=True)
    mkdir_p(case_name)
    make_solution_map_with_snapshots(archive, filename=case_name + "/map.pdf")
    # make_solution_snapshots(archive, filename=case_name + "/t={0}.pdf")
    del archive

post_process_case("lc0.01_forcing")
post_process_case("lc0.01_noforcing")
post_process_case("lc6_forcing")
post_process_case("lc6_noforcing")
