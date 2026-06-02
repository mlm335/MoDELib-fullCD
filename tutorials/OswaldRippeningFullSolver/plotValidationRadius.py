import sys, string, os
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
sys.path.insert(0, 'lib')
plt.rc('text', usetex=True)
plt.rc('font', family='Times New Roman')
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
from readNodesMod2 import *
from readFile import *
from readF import *
import re
from scipy.integrate import cumulative_trapezoid
    
def getNodeData(folderName):
    F=np.loadtxt(folderName+'/F/F_0.txt');
    print(folderName+'/F/F_0.txt has size ' + str(np.shape(F)));
    Pos=[];
    Velo=[];
    ClimbScalar=[];
    A=np.empty(np.size(F[:,0]))
    n=0;
    for runID in F[:,0]:
        nd=readNODEStxt(folderName+'/evl/evl_'+str(int(runID)))
        Pos.append(nd.nodesPos);
        Velo.append(nd.nodesV);
        ClimbScalar.append(nd.nodesClimbVScalar);
        A[n]=np.mean(nd.loopsArea);
        n=n+1;
    return F[:,1], Pos, Velo, ClimbScalar, A;
    
    
def getPeachKoehler(folderName):
    F = np.loadtxt(os.path.join(folderName, 'F', 'F_0.txt'))
    PkForce = []
    for runID in F[:, 0]:
        file_path = os.path.join(folderName, 'evl', f'ddaux_{int(runID)}.txt')
        with open(file_path, 'r') as file:
            first_number = int(file.readline().strip())
            second_number = int(file.readline().strip())
            for _ in range(first_number + 1):
                file.readline()
            magVelo = [];
            for _ in range(second_number):
                line = file.readline().strip()
                line = re.split(r'\s+', line)
                magVelo.append( np.sqrt( float(line[19])**2 + float(line[20])**2 + float(line[21])**2) )
            PkForce.append(np.mean(magVelo))
    return np.array(PkForce)
    
k_B_eV_K = 8.617333262e-5 # [eV/K]
k_B_J_K = 1.3806452e-23 #[J/K]
        
####################################################################
fileNames = ['vLoop_0','iLoop_0']
marker = ['o','v']
color = ['royalblue','firebrick']

####################################################################
#fig1, ax1 = plt.subplots(figsize=(12,8)) # Individual Velocities for Vacancy Loops
#fig2, ax2 = plt.subplots(figsize=(12,8)) # Individual Velocities for Interstitial Loops
fig3, ax3 = plt.subplots(figsize=(12,8)) # Radius vs Time for both
fig4, ax4 = plt.subplots(figsize=(12,8)) # Radius vs Time for both
####################################################################

Folder = 'MoDELib2/Validation/1000KNoRelax'
dataFolder = str(pathlib.Path().resolve()) + '/' + Folder;
print(dataFolder)

for i,k in enumerate(fileNames):
    matfilepath = dataFolder + '/' + k + '/Zr2.txt'; #directory of material file
    
    # DDD Parameters
    DoseRate = get_scalar(matfilepath,'doseRate_dpaPerSec')
    mu0_SI = get_scalar(matfilepath,'mu0_SI')
    rho_SI = get_scalar(matfilepath,'rho_SI')
    b_SI = get_scalar(matfilepath,'b_SI')
    v_dd2SI=np.sqrt(mu0_SI/rho_SI);
    t_dd2SI=b_SI/v_dd2SI;

    print("Material File Path: ", matfilepath, " for simulation ", k)
    print("Shear Wave Speed:",v_dd2SI,"Time Units DDD:",t_dd2SI)
    print("Dose Rate:",DoseRate)

    ##################### Simulation data ###############################
    Time, Pos, Velo, ClimbScalar, Area  = getNodeData(dataFolder + '/' + k)

    # Average Nodal Climb Velocity and Radius of the Loop
    nodalVClimb = []
    nodalIClimb = []
    for array in ClimbScalar:
        if array.size > 0:
            vClimb = array[:, 0]
            iCLimb = array[:, 1]
            # Filter out None or NaN values
            valid_vCLimb = [value for value in vClimb if not np.isnan(value)]
            valid_iCLimb = [value for value in iCLimb if not np.isnan(value)]
            nodalVClimb.append(np.mean(valid_vCLimb))
            nodalIClimb.append(np.mean(valid_iCLimb))
    
    v_climb_array = np.array(nodalVClimb)*v_dd2SI
    i_climb_array = np.array(nodalIClimb)*v_dd2SI
    dataLength = len(i_climb_array)
    netClimbVelo = ( v_climb_array[0:dataLength] - i_climb_array[0:dataLength] )
    timeClimb = Time[0:dataLength]*t_dd2SI
    radius = np.array ( np.sqrt(Area[0:dataLength]*b_SI*b_SI/math.pi) )
    

    ##################### Parameters for Analytical Radius ###############################
    Em_eV = get_matrix(matfilepath, 'mobileSpeciesEnergyMigration_eV') # [eV]
    Ef_eV = get_vector(matfilepath, 'mobileSpeciesEnergyFormation_eV') # [eV]
    D0_SI = get_matrix(matfilepath, 'mobileSpeciesD0_SI') # [m^2/s]
    b_SI = get_scalar(matfilepath,'b_SI') # [m]
    
    # Atomic Volume [m^3]
    ca = np.sqrt(8 / 3)
    Va = (((np.sqrt(3) * 3 / 2) * (ca) * b_SI ** 3) / 6) # * 10  # [m^3]  # ** 27  # nm^3
    r_disloc = 2*b_SI # [m]
    T = 1000 # [K]
    Dv = D0_SI[0,0] * math.exp( -Em_eV[0,0] / k_B_eV_K / T )
    Di = D0_SI[1,0] * math.exp( -Em_eV[1,0] / k_B_eV_K / T )
    cv0 = math.exp( -Ef_eV[0] / k_B_eV_K / T )
    ci0 = math.exp( -Ef_eV[1] / k_B_eV_K / T )
    mu0_SI = get_scalar(matfilepath,'mu0_SI') # [Pa/K]
    nu = get_scalar(matfilepath,'nu') # [-]

#    print('Migration Energy Vacancy [eV]:', Em_eV[0,0])
#    print('Migration Energy Interstitial [eV]:', Em_eV[1,0])
#    print('Formation Energy Vacancy [eV]:', Ef_eV[0])
#    print('Formation Energy Interstitial [eV]:', Ef_eV[1])
#    print('Isotropic Vacancy Pre-exponential Diffusion Coeff [-]', D0_SI[0,0])
#    print('Isotropic Interstitial Pre-exponential Diffusion Coeff [-]', D0_SI[1,0])
#    print('Burgers Vectore [m]',b_SI)
#    print('Volume of one atom [m^3]:',Va)
#    print('Dislocation Core Radius [m]:', r_disloc)
#    print('kb [eV/K]:', k_B_eV_K)
#    print('T [K]:', T)
#    print('Dv [m^2/s]:', Dv , 'Di [m^2/s]:', Di)
#    print('cv0 [-]:', cv0 , 'ci0 [-]:', ci0)
#    print('Loop Radius [m]:', R)
#    print('Modulus [Pa/K]:', mu0_SI)
#    print('Poisson [-]:', nu)

    ##################### Analytical Expression ###########################

    # Pk Force from Code
    PkForce = getPeachKoehler(dataFolder+'/'+k)
    pkClimb = np.array(PkForce[0:dataLength]*mu0_SI*b_SI)
    
    # Analytic peach koehler force
    fcl = [];
    for j in range(len(radius)):
        fcl.append( ( mu0_SI*b_SI*b_SI ) *  ( math.log( 8*radius[j] / r_disloc) - 1 ) / ( 4 * math.pi * (1 - nu) * radius[j] ) )
    analyticPK = np.array(fcl);

    
    b_c = np.sqrt(8/3)*b_SI/2
    # Compute without taylor expansion
    if 'iLoop' in k:  # Interstitial Loop
        Const = 2*np.pi / ( b_SI*np.log(8*radius/r_disloc) ) / 1.0
        cv_eq = cv0*np.exp( - (pkClimb*Va) / (b_SI*k_B_J_K*T) )
        ci_eq = ci0*np.exp( (pkClimb*Va) / (b_SI*k_B_J_K*T) )
        totalClimbVelocity = Const*(  Dv*cv0 - Dv*cv_eq - Di*ci0 + Di*ci_eq )
        # Plot Climb Velocity
        ax4.plot( radius*1e9, np.abs(totalClimbVelocity)*1e9, color='k')
        ax4.scatter( radius*1e9, np.abs(netClimbVelo)*1e9, color = color[i], marker=marker[i], s=90)
        # Integrate for Radius
        radius_integrated = cumulative_trapezoid(totalClimbVelocity,timeClimb,0.0)
        radius_integrated2 = np.insert(radius_integrated,0,0)
        # Plot Radius
        ax3.plot(timeClimb, (3.98984866e-08 - radius_integrated2)*1e9, color='k',linewidth=1.8)
        ax3.scatter(timeClimb[::1], radius[::1]*1e9, color = color[i], marker=marker[i], s=90)
        
    else:   # Vacancy Loop
        Const = 2*np.pi / ( b_c*np.log(8*radius/r_disloc) ) / 1.0
        cv_eq = cv0*np.exp( (pkClimb*Va) / (b_c*k_B_J_K*T) )
        ci_eq = ci0*np.exp( -(pkClimb*Va) / (b_c*k_B_J_K*T) )
        totalClimbVelocity = Const*(  Dv*cv0 - Dv*cv_eq - Di*ci0 + Di*ci_eq )
        # Plot Climb Velocity
        ax4.plot( radius*1e9, np.abs(totalClimbVelocity)*1e9, color='k')
        ax4.scatter( radius*1e9, np.abs(netClimbVelo)*1e9, color = color[i], marker=marker[i], s=90)
        # Integrate for Radius
        radius_integrated = cumulative_trapezoid(-totalClimbVelocity,timeClimb,0.0)
        radius_integrated2 = np.insert(radius_integrated,0,0)
        # Plot Radius
        ax3.plot(timeClimb,(3.98984866e-08 - radius_integrated2)*1e9, color='k',linewidth=1.8)
        ax3.scatter(timeClimb[::1], radius[::1]*1e9, color = color[i], marker=marker[i], s=90)

####################################################################
textfont = 24 #Use for x axis and y axis labels
axisfont = 20 #Use for anotations, xticks and yticks, legends
legendfont = 20

# Add the legend
legend_lines = [
    Line2D([0], [0], color='royalblue', marker=marker[0], markersize=10, linestyle='None', label=r'Vacancy $\langle c \rangle$ Loop: Thermal'),
    Line2D([0], [0], color='firebrick', marker=marker[1], markersize=10, linestyle='None', label=r'Interstitial $\langle a \rangle$ Loop: Thermal'),
    Line2D([0], [0], color='k', lw=2, label=r'Analytical Radius: Thermal')
]
legend1_ax3 = ax3.legend(handles=legend_lines, fontsize=legendfont, loc='upper right')
ax3.add_artist(legend1_ax3)
ax3.set(xlabel=r'Time [sec]', ylabel=r'Loop Radius [nm]');
ax3.grid(True, linestyle='dotted')
for item in ([ax3.xaxis.label, ax3.yaxis.label]): item.set_fontsize(textfont)
ax3.xaxis.offsetText.set_fontsize(legendfont)
ax3.yaxis.offsetText.set_fontsize(legendfont)
ax3.tick_params(axis = 'both', which = 'major', direction='in', bottom=True, left=True, top=True, right=True, labelsize = axisfont)
ax3.tick_params(axis = 'both', which = 'minor', direction='in', bottom=True, left=True, top=True, right=True, labelsize = axisfont)
ax3.set_xlim(0)
ax3.set_ylim(5)


ax4.invert_xaxis()
ax4.set(xlabel=r'Loop Radius [nm]', ylabel=r'Climb Velocity [nm/s]');
ax4.grid(True, linestyle='dotted')
for item in ([ax4.xaxis.label, ax4.yaxis.label]): item.set_fontsize(textfont)
ax4.xaxis.offsetText.set_fontsize(legendfont)
ax4.yaxis.offsetText.set_fontsize(legendfont)
ax4.tick_params(axis = 'both', which = 'major', direction='in', bottom=True, left=True, top=True, right=True, labelsize = axisfont)
ax4.tick_params(axis = 'both', which = 'minor', direction='in', bottom=True, left=True, top=True, right=True, labelsize = axisfont)
ax4.set_xlim(40,1.5)
ax4.set_ylim(-0.005,0.06)

fig3.savefig('/Users/matthewmaron/Documents/MoDELib2-CD/Results/DiscreteGrowth/Figures/ValidationRadiusNew.pdf', format='pdf', dpi=1500)
#fig4.savefig('/Users/matthewmaron/Documents/MoDELib2-CD/Results/DiscreteGrowth/Figures/ValidationVelocity.pdf', format='pdf', dpi=1500)

plt.show()


