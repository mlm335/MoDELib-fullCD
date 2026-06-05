import numpy as np
import sys, os
from matplotlib import pyplot as plt
print(os.environ['PATH'])
sys.path.append('../../python/')
sys.path.append('../../lib/')

from readF import *
from readFile import *
from readEVL import *
import pathlib

# -- Experimental Data Carpenter Rogerson -- #
EXPIa=np.array([[0.1245,1.3249],[0.3034,1.1033],[0.6301,1.163],[0.8557,1.0831],
    [1.3613,1.4861],[1.5946,1.3451],[1.7191,1.4055],[2.5514,1.3048],[2.7381,1.2645],
    [3.7105,1.6474],[4.4417,2.0101],[4.6283,2.5945],[5.7174,2.9572],
    [6.0596,2.9773],[6.1141,3.1184],[6.6197,3.4207],[7.2887,4.3073],[8.1677,5.1537]])  #IODIDE PURITY a # x dpa y strain
EXPIc=np.array([[ 0.1400,-0.1461],[0.3034,-0.4081],[0.6379, -0.4081],[1.2368,-0.7103],
    [1.5013,-0.8111],[1.6180,-0.6902],[2.1936,-0.6902],[2.6059,-1.3149],
    [3.1037,-1.1537],[3.5860,-1.2544],[4.2083,-1.0126]]) #IODIDE PURITY c # x dpa y strain
EXPZa=np.array([[0.1245,0.5189],[0.2723,0.7406],[0.4745,0.8816],[0.6845,1.1234],[1.0346,1.4660],[1.2757,1.5063],
    [1.5013,1.2645],[1.6024,1.0428],[1.9758,1.4660],[2.3492,1.2645],[3.1271,1.7884],
    [3.3293,1.8690],[3.6482,2.1310],[3.9594,2.1713],[4.4261,2.3325],[4.6283,2.1713],
    [5.6940,2.4332],[5.9507,2.9370], [6.4719,2.8766],[6.6430,3.1385],[6.7753,3.1788]])  #Zone Refined a # x dpa y strain units e-4
EXPZc=np.array([[0.1478,-0.4282],[0.2800,-0.2267],[0.6145,-0.7909],[0.8634,-0.7506],[1.3613,-0.7506],
    [1.4935,-0.6297],[1.8124,-0.7708],[2.1625,-0.8917],[2.3647,-0.7103],[3.4693,-1.0126]]) #Zone Refined c # x dpa y strain units e-4

# ------------------------------------------------------------------------ #
#                          Main Inputs
# ------------------------------------------------------------------------ #
folder = str(pathlib.Path().resolve()) # Adjust this path as needed to directory with simulation folders
textfont = 18
axisfont = 16

# ------------------------------------------------------------------------ #
#                          Main Processing Loop
# ------------------------------------------------------------------------ #
print(f"Processing {folder}")
matfile = os.path.join(folder, 'inputFiles', 'Zr4_Fitted.txt')
F, Flabels = readFfile(folder)
t = getFarray(F, Flabels, 'time [b/cs]')
DoseRate = get_scalar(matfile, 'doseRate_dpaPerSec')
mu0_SI = get_scalar(matfile, 'mu0_SI')
rho_SI = get_scalar(matfile, 'rho_SI')
b_SI = get_scalar(matfile, 'b_SI')
v_dd2SI = np.sqrt(mu0_SI / rho_SI)
t_dd2SI = b_SI / v_dd2SI
dpa = t * t_dd2SI * DoseRate  # <-- missing before

# Plotting strain components and swelling
b11 = getFarray(F, Flabels, 'betaP_11')
b22 = getFarray(F, Flabels, 'betaP_22')
b33 = getFarray(F, Flabels, 'betaP_33')
b12 = getFarray(F, Flabels, 'betaP_12')
b13 = getFarray(F, Flabels, 'betaP_13')
b23 = getFarray(F, Flabels, 'betaP_23')

swelling_trace_over_3 = []
for i in range(len(t)):
    betaP = np.array([
        [b11[i], b12[i], b13[i]],
        [b12[i], b22[i], b23[i]],
        [b13[i], b23[i], b33[i]],
    ])
    trace = np.trace(betaP)
    swelling_trace_over_3.append(trace / 3)
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(dpa+0.5, b11*100, label=r'$\beta_{11}$', color='red')
ax.plot(dpa+0.5, b22*100, label=r'$\beta_{22}$', color='green')
ax.plot(dpa+0.5, b33*100, label=r'$\beta_{33}$', color='blue')
ax.plot(dpa+0.5, swelling_trace_over_3, label=r'Swelling ($\mathrm{tr}(\beta^P)/3$)', color='gray', linestyle='-')


ax.set_xlim(0,9.2)
ax.set_ylim(-0.08,0.06)

# Experimental Data
ax.scatter(EXPIa[:,0],EXPIa[:,1]*10**-2,color="firebrick",marker="s",s=80,label=r'Experiment: Iodide Purity');
ax.scatter(EXPZa[:,0],EXPZa[:,1]*10**-2,color="orange",marker="s",s=80,label=r'Experiment: Zone Refined Purity');
ax.scatter(EXPIc[:,0],EXPIc[:,1]*10**-2,color="firebrick",marker="s",s=80);
ax.scatter(EXPZc[:,0],EXPZc[:,1]*10**-2,color="orange",marker="s",s=80);


ax.set_xlabel('Irradiation Dose [dpa]', fontsize=textfont)
ax.set_ylabel('Strain [%]', fontsize=textfont)
ax.legend(fontsize=axisfont)
ax.grid(True, linestyle=':')
ax.set_xlim(left=0)
ax.axhline(0, color='k', linewidth=0.5)
ax.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig.tight_layout()



# Plot dislocation density
glissile_density = getFarray(F, Flabels, 'glissile density [m^-2]')
sessile_density  = getFarray(F, Flabels, 'sessile density [m^-2]')
fig2, ax2 = plt.subplots(figsize=(10, 6))
ax2.plot(dpa, glissile_density+sessile_density, label='Dislocation Density', color='k', linewidth=2)
ax2.set_xlabel('Irradiation Dose [dpa]', fontsize=textfont)
ax2.set_ylabel('Dislocation Density [m$^{-2}$]', fontsize=textfont)
ax2.legend(fontsize=axisfont)
ax2.grid(True, linestyle=':')
ax2.set_xlim(left=0)
ax2.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig2.tight_layout()

plt.show()
# Save figures
save_path = os.path.join(folder, 'StrainComponentsWithSwelling.pdf')
fig.savefig(save_path)
save_path2 = os.path.join(folder, 'DislocationDensities_vs_DPA.pdf')
fig2.savefig(save_path2)   
plt.close(fig)
print(f"Saved {save_path}")
plt.close(fig2)
print(f"Saved {save_path2}")

