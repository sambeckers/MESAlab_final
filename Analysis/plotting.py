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

# HR Diagram
def plot_hr_diagram():
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

l = mr.MesaLogDir('/Users/sam/Documents/GitHub/MESAlab_final/2M_prems_to_WD/LOGS')
prof_nums = l.profile_numbers

# Core evolution
def plot_core_evolution():

    central_Ts = []
    central_Rhos = []

    for prof in prof_nums:
        n = l.profile_data(profile_number=prof)
        central_Ts.append(n.logT[-1])
        central_Rhos.append(n.logRho[-1])

    plt.figure(dpi=300)
    plt.plot(central_Ts, central_Rhos)
    plt.show()

def plot_convective_preMS():
    # Define the parameters
    # X = 0.7
    # Y = 0.28
    # Z = 0.02

    # # Mean molecular weights
    # mu = 1 / ((2 * X) + (0.75 * Y) + (0.5 * Z))
    # mu_e = 2 / (1 + X)

    # # Temperature range
    # T = np.linspace(6, 8, 258)

    # # Ideal gas density
    # rho_ideal = (10**(T) / (3.2e7 * (mu**(-1/3))))**3
    
    # plt.plot(T, rho_ideal)
    # plt.show()

    n = l.profile_data(profile_number=50)
    grad_R = n.gradr
    grad_A = n.grada
    radius = n.radius

    plt.figure(dpi=300)
    plt.plot(radius, grad_R, label='Radiative Gradient')
    plt.plot(radius, grad_A, label='Convective Gradient')
    plt.show()

def plot_convective_MS():
    pass

if __name__ == "__main__":
    # plot_hr_diagram()
    # plot_core_evolution()
    # plot_convective_preMS()
    # plot_convective_MS()