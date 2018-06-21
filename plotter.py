import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here,'./python/')))
from matplotlib.backends.backend_pgf import FigureCanvasPgf
from matplotlib.backend_bases import register_backend
register_backend('pdf', FigureCanvasPgf)
import planckStyle as s
from pylab import *


from getdist import plots, MCSamples
import getdist

print('Version: ',getdist.__version__)

import GetDistPlots

import planckStyle

g = planckStyle.getSinglePlotter(chain_dir = './chains', ratio=1)


roots = ['JLA+BAO_4bins_binned','JLA+BAO_4bins_binned_originalbins']
params = ['omegam','binw1','binw2','binw3','binw4']
colors = ['#8E001C','#FFB300']
g.settings.solid_contour_palefactor = 0.8
g.triangle_plot(roots, params, filled=[True,True,False], contour_colors=colors, legend_colors=colors, legend_labels=[r'new bins','old bins','tau'], legend_loc='upper right',tight_layout=True)


g.export('results_plots/triangle.pdf')
