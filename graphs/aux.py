import numpy as np
import matplotlib as mpl
# Define a backend
# MUST be invoked before pyplot import
mpl.use("pgf")
import matplotlib.pyplot as plt

# For R-style plotting
plt.style.use('ggplot')


def set_plt_params(
        relative_fig_width=1.0, landscape=True, page_width=307.3, rescale_h=1):

    fig_width_pt = page_width * relative_fig_width
    inches_per_pt = 1.0 / 72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches

    if landscape:
        fig_height = fig_width * golden_mean  # height in inches
    else:
        fig_height = fig_width / golden_mean  # height in inches

    fig_height = fig_height * rescale_h
    fig_size = [fig_width, fig_height]
    params = {
        'font.family': 'serif',
        'axes.labelsize': 7,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'axes.labelcolor': 'black',
        'ytick.color': 'black',
        'xtick.color': 'black',
        'legend.handlelength': 4,
        'legend.fontsize': 7,
        # 'lines.markersize': 3,
        # 'xtick.labelsize': 7,
        # 'ytick.labelsize': 7,
        'text.usetex': True,
        'text.latex.unicode': True,
        'figure.figsize': fig_size,
        'pgf.texsystem': "xelatex",
        'pgf.rcfonts': False,
    }

    plt.rcParams.update(params)

plt.ioff()
