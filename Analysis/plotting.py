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
star_age = []

# Core evolution
def plot_core_evolution():

    central_Ts = []
    central_Rhos = []
    star_age = []

    for prof in prof_nums:
        n = l.profile_data(profile_number=prof)
        central_Ts.append(n.logT[-1])
        central_Rhos.append(n.logRho[-1])
        star_age.append(n.star_age)

    star_age = np.array(star_age)

    # Solar abundances
    X = 0.7
    Y = 0.28
    Z = 0.02

    #mean molecular weights
    mu = 1/((2*X)+(0.75*Y)+(0.5*Z))
    mu_e = 2/(1+X)

    # Constants
    a = 7.56578e-15 # erg cm^-3 K^-4
    R = 8.31447e7 # erg g^-1 K^-1
    K_NR = 1.0036e13 #cgs
    K_ER = 1.2435e15 #cgs

    logT = np.linspace(np.min(central_Ts)-1, np.max(central_Ts)+1,258)

    # Radiation - ideal gas density boundary
    rho_radid = (mu * a * (10**logT)**3) / (3*R)

    # Ideal gas - non-degenerate NR density boundary
    rho_idER = ((mu_e**(5/3) * R * (10**logT)) / (K_NR * mu))**(3/2)

    # non-degenerate NR - degenerate ER density boundary
    rho_nrer = (K_ER/K_NR)**3 * mu_e

    # Theoretical evolution
    G = 6.6743 * 1e-8 
    M = 1.9884 * 1e33
    rho_ideal_theory = ((10**(logT)*R)/(G*mu*(M**(2/3))))**3
    rho_NR_theory = (G/K_NR)**3 * mu_e**5 * M**2


    plt.figure(dpi=300, figsize=(10, 8))
    plt.plot(central_Ts, central_Rhos, c='grey', alpha=0.5, zorder=1)
    sc = plt.scatter(central_Ts, central_Rhos, s=20, c=star_age/1e9, cmap='viridis', zorder=2)
    plt.scatter(np.log10(1.5e7), np.log10(150), s=100, c='red', marker='*', label='Sun', zorder=3)
    plt.fill_between(logT, -8, np.log10(rho_radid), color='purple', alpha=0.4, zorder=0)
    plt.plot(logT, np.log10(rho_radid), c='k', linestyle = '--', zorder=0)
    plt.fill_between(logT, np.log10(rho_radid), np.log10(rho_idER), color='purple', alpha=0.2, zorder=0)
    plt.plot(logT, np.log10(rho_idER),  c='k', linestyle = '--', zorder=0)
    plt.fill_between(logT, np.log10(rho_idER), np.log10(rho_nrer), color='purple', alpha=0.1, zorder=0)
    plt.axhline(y=np.log10(rho_nrer), color='k', linestyle='--', zorder=0)
    plt.fill_between(logT, np.log10(rho_nrer), 8.5, color='purple', alpha=0.05, zorder=0)
    plt.plot(logT, np.log10(rho_ideal_theory), c='r', linestyle = '--', zorder=0, label='Theoretical evolution')
    plt.axhline(y=np.log10(rho_NR_theory), c='r', linestyle = '--', zorder=0)
    plt.xlim(5, 9)
    plt.ylim(-7, 8)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.xlabel('$\log_{10}{T_{\mathrm{c}}}$ [K]', fontsize=16)
    plt.ylabel('$\log_{10}{\\rho_{\mathrm{c}}}$ [g cm$^{-3}$]', fontsize=16)
    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Star Age (Gyr)', fontsize=16)
    plt.legend(fontsize=16, loc='lower right')
    plt.savefig('/Users/sam/Documents/GitHub/MESAlab_final/Analysis/core_evolution.svg')
    plt.show()

def plot_convective_preMS():
    n = l.profile_data(profile_number=7)
    grad_R = n.gradr
    grad_A = n.grada
    radius = n.logR
    print('Age of the star (MS):', n.star_age/10**9, 'Gyr')
    index = np.argmax(grad_R > grad_A)
    radius_value = radius[index]

    print(f"First radius where gradR > gradA: 10^{radius_value:.3f} cm")
    plt.figure(dpi=300)
    plt.plot(radius, grad_R, c='slateblue', ls = '--', lw=2, label='Radiative Gradient')
    plt.plot(radius, grad_A, c='darkblue', ls = '--', lw=2, label='Adiabatic Gradient')
    # plt.axvline(x=np.log10(radius_value), color='blue', ls = '--', label='Convective Boundary') 
    ylims = plt.ylim()
    plt.fill_betweenx(np.linspace(*ylims, 500), np.min(radius), np.log10(radius_value), color='grey', alpha=0.2, label='Convective Zone')
    plt.yscale('log')
    plt.ylim(np.min(grad_R),np.max(grad_R)+0.05)
    plt.xlim(np.min(radius))
    plt.xlabel(r'$\log{R} [R_{\odot}]$', fontsize=16)
    plt.ylabel(r'$\nabla$', fontsize=16)
    plt.legend()
    plt.savefig('/Users/sam/Documents/GitHub/MESAlab_final/Analysis/convective_preMS.svg')
    plt.show()

def plot_convective_MS():
    n = l.profile_data(profile_number=8)
    grad_R = n.gradr
    grad_A = n.grada
    radius = n.logR
    print('Age of the star (MS):', n.star_age/10**9, 'Gyr')
    index = np.argmax(grad_R > grad_A)
    radius_value = radius[index]

    print(f"First radius where gradR > gradA: 10^{radius_value:.3f} cm")
    plt.figure(dpi=300)
    plt.plot(radius, grad_R, c='slateblue', ls = '--', lw=2, label='Radiative Gradient')
    plt.plot(radius, grad_A, c='darkblue', ls = '--', lw=2, label='Adiabatic Gradient')
    # plt.axvline(x=np.log10(radius_value), color='blue', ls = '--', label='Convective Boundary') 
    ylims = plt.ylim()
    plt.fill_betweenx(np.linspace(*ylims, 500), np.min(radius), np.log10(radius_value), color='grey', alpha=0.2, label='Convective Zone')
    plt.yscale('log')
    plt.ylim(np.min(grad_R),np.max(grad_R)+0.05)
    plt.xlim(np.min(radius))
    plt.xlabel(r'$\log{R} [R_{\odot}]$', fontsize=16)
    plt.ylabel(r'$\nabla$', fontsize=16)
    plt.legend()
    plt.savefig('/Users/sam/Documents/GitHub/MESAlab_final/Analysis/convective_MS.svg')
    plt.show()

if __name__ == "__main__":
    plot_hr_diagram()
    plot_core_evolution()
    plot_convective_preMS()
    plot_convective_MS()