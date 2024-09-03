import argparse
import h5py as h5
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

cmap = cm.magma

parser = argparse.ArgumentParser(prog='preview.py')
parser.add_argument('filename')
args = parser.parse_args()

fig, axs = plt.subplots(nrows=2, height_ratios=[0.1, 1])

with h5.File(args.filename, 'r') as dat:
    rho = dat['rho']
    rho_ux = dat['rho_ux']
    rho_uy = dat['rho_uy']

    init, nstep = 0, rho.shape[0]
    slider = Slider(axs[0], 'step', 0, nstep - 1,
        valinit=init, valstep=list(range(nstep)))

    x = np.linspace(0, 1, rho.shape[1])
    y = np.linspace(0, 1, rho.shape[2])

    def update(val):
        Q.set_array(rho[val])
        Q.autoscale()

        cbar.update_normal(Q)
        fig.canvas.draw_idle()

    Q = axs[1].pcolormesh(x, y, rho[init], cmap=cmap)
    cbar = fig.colorbar(Q, ax=axs[1])
    slider.on_changed(update)

    axs[1].set_aspect(1)
    plt.show()

# vim: set ft=python:
