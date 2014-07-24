import numpy as np
from matplotlib.pyplot import plot, show, legend, close
import matplotlib.pyplot as plt
close("all")
import seaborn as sns
#sns.set(style="whitegrid")
sns.set(style="ticks")

# Allow text to be edited in Illustrator
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

res = np.load("SphereResults.npy")
#res = res['arr_0']
radii=[85,88,92,100]
max_r = np.max(radii)
bound = 82 # min(radii)
step = 4

rdm = False
mags = True

resSimbio = np.load("SimBioRDM_MAG.npy")

if rdm:
    plot(res[range(0,bound,step),0], res[range(0,bound,step),1], linewidth=2, label="GetDP")
    plot(res[range(0,bound,step),0], res[range(0,bound,step),3], linewidth=2, label="OpenMEEG")
    plot(res[range(0,bound,step),0], resSimbio[range(0,bound,step),0][::-1], linewidth=2, label="NeuroFEM/SimBio")
    plt.xlabel('Dipole position (from 0 in Z dir, max radius = %3.0f)' % max_r, fontsize=14)
    plt.ylabel('Relative difference measure (RDM)', fontsize=14)
    for r in radii:
        plt.axvline(x=r,ymin=0,ymax=1,color="k",linestyle="--")
    legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=2, mode="expand", borderaxespad=0.)
    plt.axhline(y=0,xmin=0,xmax=100, color="r",linestyle="--")
    show()
    plt.axis([10, 105, -0.02, 0.15])
    ax = plt.gca()
    ax.set_autoscale_on(False)
    plt.savefig("RDMs.pdf")

#close("all")

if mags:
    # Plot MAGs
    #plot(res[range(0,82,6),0], res[range(0,82,6),1], linewidth=2, label="GetDP")
    plot(res[range(0,bound,step),0], res[range(0,bound,step),2], linewidth=2, label="GetDP")
    plot(res[range(0,bound,step),0], res[range(0,bound,step),4], linewidth=2, label="OpenMEEG")
    plot(res[range(0,bound,step),0], resSimbio[range(0,bound,step),1][::-1], linewidth=2, label="NeuroFEM/SimBio")
    plt.xlabel('Dipole position (from 0 in Z, max radius = %3.0f)' % max_r, fontsize=14)
    plt.ylabel('Magnification errors (MAG)', fontsize=14)
    for r in radii:
        plt.axvline(x=r,ymin=0,ymax=1,color="k",linestyle="--")

    plt.axhline(y=1,xmin=0,xmax=100, color="r",linestyle="--")
    legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=2, mode="expand", borderaxespad=0.)
    show()
    plt.axis([10, 105, 0.5, 1.2])
    ax = plt.gca()
    ax.set_autoscale_on(False)
    plt.savefig("MAGs.pdf")

#close("all")
