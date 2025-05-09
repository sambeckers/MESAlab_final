import mesa_reader as mr
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman",
    "font.sans-serif": "helvetica"
})
history = mr.MesaData('/Users/sam/Documents/GitHub/MESAlab_final/2M_prems_to_WD/LOGS/history.data')

plt.figure(dpi=300, figsize=(10, 6))
plt.scatter(history.log_Teff, history.log_L, s=5, c=history.star_age/10**9, cmap='jet')
plt.tick_params(axis='both', which='major', labelsize=14)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Star Age (Gyr)', fontsize=14)
plt.xlabel('$\log{T_{\mathrm{eff}}}$ [K]', fontsize=14)
plt.ylabel('$\log{L} [L_{\odot}$]', fontsize=14)
plt.gca().invert_xaxis()
plt.grid(linestyle='--', alpha=0.7)
plt.savefig('/Users/sam/Documents/GitHub/MESAlab_final/Analysis/HR.svg')
plt.show()
