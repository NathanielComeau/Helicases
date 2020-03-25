# Try to plot everything in one big happy plot!

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
import scipy.stats

fig = plt.figure(figsize=(6,10))
ax1 = fig.add_axes([0.0, 0.8, 1.0, 0.2])
ax1.set_xticks([])
ax1.set_yticks([])
ax2 = fig.add_axes([0.0, 0.6, 1.0, 0.2])
ax2.set_xticks([])
ax2.set_yticks([])
ax3 = fig.add_axes([0.0, 0.4, 1.0, 0.2])
ax3.set_xticks([])
ax3.set_yticks([])
ax4 = fig.add_axes([0.0, 0.2, 1.0, 0.2])
ax4.set_xticks([])
ax4.set_yticks([])
ax5 = fig.add_axes([0.0, 0.0, 1.0, 0.2])
ax5.set_xticks([])
ax5.set_yticks([])

def myplot(x, y, w, s, b):
    heatmap, xedges, yedges, k = scipy.stats.binned_statistic_2d(x, y, w, bins=b)
    med = np.median(heatmap)
    heatmap = np.nan_to_num(heatmap, nan=med)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent

axs = [ax1, ax2, ax3, ax4, ax5]

# We'll need the max and min of all our data to set a common colour scale
_min, _max, _median = 0, 0, 0
# Use min and max from first frame as colour bar for all frames
with open('bin1_qa.txt') as f:
    lines = f.readlines()
average_qas = [float(x.rstrip('\n')) for x in lines]
_min = min(_min, min(average_qas))
_max = max(_max, max(average_qas))

for ax, i in zip(axs, range(1,6)):
    # Load up three parameters
    with open('bin'+str(i)+'_x_coord.txt') as f:
        lines = f.readlines()
    x_lines = [int(x.rstrip('\n')) for x in lines]
    with open('bin'+str(i)+'_y_coord.txt') as f:
        lines = f.readlines()
    y_lines = [int(x.rstrip('\n')) for x in lines]
    with open('bin'+str(i)+'_qa.txt') as f:
        lines = f.readlines()
    average_qas = [float(x.rstrip('\n')) for x in lines]

    img, extent = myplot(x_lines, y_lines, average_qas, 0, 100)
    res = ax.imshow(img, extent=extent, origin='lower', cmap=cm.jet, vmin=_min, vmax=_max)
    #ax.set_title(r"Averaging and binning  $\sigma = %d, \beta$ = %d" % (s,b))
    #ax.set_title(r"Averaging and binning, $\beta$ = %d" % (s))
    #fig.colorbar(res, ax=ax)

fig.suptitle('2D average binning in x and y of quality scores', fontsize=16)

#plt.tight_layout()
#plt.savefig('Entire_Human_Genome_x_y_qa.png')
plt.show()
