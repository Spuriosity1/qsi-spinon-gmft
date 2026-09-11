#!/bin/env python3
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse


tex_fonts_ss = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "sans-serif",
    "text.latex.preamble": r"\usepackage{sfmath}\renewcommand{\familydefault}{\sfdefault}",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}


tex_fonts_serif = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}


plt.rcParams.update(tex_fonts_serif)


ap = argparse.ArgumentParser('plot_hhl')
ap.add_argument('file', help="hdf5 fiel contianing HHL data")
ap.add_argument('dataset', 
                choices=('Spm', 'Spp', 'Smag'), default='Spm',
                help="hdf5 fiel contianing HHL data")
ap.add_argument('--cmap',
                choices=list(mpl.colormaps),
                default='plasma'
                )

a = ap.parse_args()

f = h5py.File(a.file, 'r')

h = f['h']
l = f['l']
data = f[a.dataset]

TITLES = dict(
        Spm=r'$\langle S^+(q) S^-(-q) \rangle$',
        Spp=r'$\langle S^+(q) S^+(-q) \rangle$',
        Smag=r'$\langle \mathbf{m}(q)\mathbf{m}(-q) \rangle$'
        )

fig, ax = plt.subplots()
mesh = ax.pcolormesh(h, l, data, shading='nearest', cmap=a.cmap)
fig.colorbar(mesh)
ax.set_title(TITLES[a.dataset])
ax.set_xlabel(r'$(h,h,\#)$')
ax.set_ylabel(r'$(\#,\#,l)$')
plt.show()
