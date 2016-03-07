import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
mpl.rc('font', size=20.)
mpl.rc('font', family='serif')
mpl.rc('text', usetex=True)
import pyGadget

plt.ioff()
sim = pyGadget.sim.Simulation('stampede/vanilla')
sim.refine_by_mass(False)
sim.set_coordinate_system('physical')
snap = sim.load_snapshot(1900)

scales = ['140 kpc (comoving)', '1 kpc (physical)', '10 pc (physical)', '100 pc (physical)']
ratio = [.1788, .1, None, .1]
zoom = ['right', 'down', None, 'left']
clims = [(-2.5,1.5), (-2.,2.), (1.5,7.5), (-0.5,5.)]
ticks = [(-2,-1,0,1), (-1,0,1), (2,3,4,5,6,7), (0,1,2,3,4)]
cpad = [-17, -17, -15, -16]
clabel = [False, True, False, True]
bbox_props = dict(boxstyle="round", fc="k", ec="k", alpha=0.5)
zc = 'w'
zls = '--'
zlw = 1.5

count = 0
for theta in np.linspace(0,2*np.pi, 250):
    view = [('y',theta)]
    imlist = []
    for scale in ['5592pc', '1000pc', '10pc', '100pc']:
        imlist.append(pyGadget.visualize.project(snap, 'ndensity', 
                                                 scale, view, centering='avg'))

    fig = plt.figure(1, (12., 12.), dpi=600)
    grid = ImageGrid(fig, 111, # similar to subplot(111)
                    nrows_ncols = (2, 2), # creates 2x2 grid of axes
                    axes_pad=0.0, # pad between axes in inch.
                    cbar_mode = 'each', cbar_size='7%', cbar_pad=0.
                    )

    for i in range(4):
        x = imlist[i][0]
        y = imlist[i][1]
        im = imlist[i][2]
        ax = grid[i]
        img = ax.imshow(im, cmap=plt.cm.bone, origin='lower')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        img.set_clim(clims[i])

        cb = plt.colorbar(img, cax=grid.cbar_axes[i])
        cb.set_ticks(ticks[i])
        cb.ax.tick_params(left='on', pad=cpad[i],
                          labelsize=15, labelcolor='k', labelleft='on', labelright='off')
        if clabel[i]: cb.set_label('Log Number Density [cm$^{-3}$]')

        ax.text(0.5, 0.025, scales[i], color='w', ha='center', va='bottom', size=12, 
                transform=grid[i].transAxes, bbox=bbox_props)

        if ratio[i]:
            axmin, axmax = ax.get_xlim()
            axlength = axmax - axmin
            mid = axlength/2
            s = ratio[i] * axlength
            s00 = [mid - s/2, mid - s/2]
            s01 = [mid - s/2, mid + s/2]
            s11 = [mid + s/2, mid + s/2]
            ax.add_line(plt.Line2D(s00, s01, c=zc, lw=zlw))
            ax.add_line(plt.Line2D(s11, s01, c=zc, lw=zlw))
            ax.add_line(plt.Line2D(s01, s00, c=zc, lw=zlw))
            ax.add_line(plt.Line2D(s01, s11, c=zc, lw=zlw))
            if zoom[i] == 'right':
                ax.add_line(plt.Line2D([mid+s/2, axmax], [mid+s/2, axmax], c=zc, lw=zlw, ls=zls))
                ax.add_line(plt.Line2D([mid+s/2, axmax], [mid-s/2, axmin], c=zc, lw=zlw, ls=zls))
            elif zoom[i] == 'down':
                ax.add_line(plt.Line2D([mid-s/2, axmin], [mid-s/2, axmin], c=zc, lw=zlw, ls=zls))
                ax.add_line(plt.Line2D([mid+s/2, axmax], [mid-s/2, axmin], c=zc, lw=zlw, ls=zls))
            elif zoom[i] == 'left':
                ax.add_line(plt.Line2D([mid-s/2, axmin], [mid+s/2, axmax], c=zc, lw=zlw, ls=zls))
                ax.add_line(plt.Line2D([mid-s/2, axmin], [mid-s/2, axmin], c=zc, lw=zlw, ls=zls))

    fig.savefig(sim.plotpath+'/'+sim.name+'/rotations/box/{:0>4}-zoom.png'.format(count), 
                bbox_inches='tight', dpi=100)
    count += 1
