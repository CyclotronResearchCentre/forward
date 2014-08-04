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
resSimbio = np.load("SimBioRDM_MAG.npy")

radii=[85,88,92,100]
max_r = np.max(radii)
bound = 85 # min(radii)
step = 4

depth = res[range(0,bound,step),0]
rdm_getdp = res[range(0,bound,step),1]
rdm_openmeeg = res[range(0,bound,step),3]
rdm_simbio = resSimbio[range(0,bound,step),0][::-1]

mag_getdp = res[range(0,bound,step),2]
mag_openmeeg = res[range(0,bound,step),4]
mag_simbio = resSimbio[range(0,bound,step),1][::-1]

print("Mean")
print("RDM, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.mean(rdm_getdp), np.mean(rdm_openmeeg), np.mean(rdm_simbio)))
print("MAGs, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.mean(mag_getdp), np.mean(mag_openmeeg), np.mean(mag_simbio)))

print("Standard Deviation")
print("RDM, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.std(rdm_getdp), np.std(rdm_openmeeg), np.std(rdm_simbio)))
print("MAGs, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.std(mag_getdp), np.std(mag_openmeeg), np.std(mag_simbio)))

#### 
print("GetDP")
print("RDM, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.mean(rdm_getdp), np.mean(rdm_openmeeg), np.mean(rdm_simbio)))
print("MAGs, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.mean(mag_getdp), np.mean(mag_openmeeg), np.mean(mag_simbio)))

print("Standard Deviation")
print("RDM, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.std(rdm_getdp), np.std(rdm_openmeeg), np.std(rdm_simbio)))
print("MAGs, GetDP: %3.3f, OpenMEEG: %3.3f, SimBio: %3.3f" % (np.std(mag_getdp), np.std(mag_openmeeg), np.std(mag_simbio)))



# This is for copy/pasting into LaTeX
print("Method & RDM & MAG\\\\")
print("\hline")
print("GetDP & %3.3f (%3.3f) & %3.3f (%3.3f)\\\\" % (np.mean(rdm_getdp), np.std(rdm_getdp), np.mean(mag_getdp), np.std(mag_getdp)))
print("\hline")
print("OpenMEEG & %3.3f (%3.3f) & %3.3f (%3.3f)\\\\" % (np.mean(rdm_openmeeg), np.std(rdm_openmeeg), np.mean(mag_openmeeg), np.std(mag_openmeeg)))
print("\hline")
print("SimBio & %3.3f (%3.3f) & %3.3f (%3.3f)\\\\" % (np.mean(rdm_simbio), np.std(rdm_simbio), np.mean(mag_simbio), np.std(mag_simbio)))