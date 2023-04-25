"""
Rayleigh-Lamb Equation Solver
Update 2023/04/25
written by: Armin Huber, armin.huber@dlr.de
-------------------------------------------------------------------------------
Calculate Lamb wave dispersion diagrams for a free isotropic plate by
evaluating the Rayleigh-Lamb equations. The Rayleigh-Lamb equation amplitudes
are scaled to a grayscale. Modal solutions are represented by minima, i.e., by
dark shading.
"""

#%% settings
import numpy as np
import matplotlib.pyplot as plt

MaterialName = 'aluminum'
MaterialLongitudinalVelocity = 6320 # (m/s)
MaterialTransverseVelocity = 3130 # (m/s)

Thickness = 1 # plate thickness (mm)

PhaseVelocityLimit = 20 # (m/ms)
FrequencyLimit = 10000 # (kHz)
PhaseVelocitySteps = 1000
FrequencySteps = 1000

#%% computation
Half = Thickness/2e3
Frequency = np.arange(0,FrequencyLimit+FrequencyLimit/FrequencySteps,FrequencyLimit/FrequencySteps);Frequency[0] = 1 # generate a range of frequencies
PhaseVelocity = np.arange(0,PhaseVelocityLimit*1e3+PhaseVelocityLimit*1e3/PhaseVelocitySteps,PhaseVelocityLimit*1e3/PhaseVelocitySteps);PhaseVelocity[0] = 1 # generate a range of phase velocities
AngularFrequency = 2*np.pi*Frequency*1e3
k2 = (AngularFrequency/np.matrix.transpose(np.tile(PhaseVelocity,(AngularFrequency.size,1))))**2 # wavenumber^2 of Lamb waves
kL2 = np.tile((AngularFrequency/MaterialLongitudinalVelocity)**2,(PhaseVelocity.size,1)) # wavenumber^2 of longitudinal bulk waves
kT2 = np.tile((AngularFrequency/MaterialTransverseVelocity)**2,(PhaseVelocity.size,1)) # wavenumber^2 of transverse bulk waves
x = np.lib.scimath.sqrt(kL2-k2) # out-of-plane wavenumber component of longitudinal bulk waves
y = np.lib.scimath.sqrt(kT2-k2) # out-of-plane wavenumber component of transverse bulk waves
a1 = (y**2-k2)**2
a2 = 4*k2*x
a3 = np.tan(x*Half)
a4 = np.tan(y*Half)
S = np.abs(a1/a3/y+a2/a4) # Rayleigh-Lamb equation for symmetric modes (absolute value)
A = np.abs(a1*a3/y+a2*a4) # Rayleigh-Lamb equation for antisymmetric modes (absolute value)

# %% dispersion diagrams
plt.figure(figsize=[16,10])
plt.imshow(20*np.log10(S),extent=[0,FrequencyLimit,0,PhaseVelocityLimit],cmap='gray',origin='lower',aspect='auto')
plt.title('Symmetric mode dispersion diagram of '+str(Thickness)+' mm '+MaterialName,fontsize=20) # Rayleigh-Lamb equation amplitude in dB
plt.xlabel('Frequency (kHz)',fontsize=20)
plt.ylabel('Phase velocity (m/ms)',fontsize=20)
plt.tick_params(labelsize=15)

plt.figure(figsize=[16,10])
plt.imshow(20*np.log10(A),extent=[0,FrequencyLimit,0,PhaseVelocityLimit],cmap='gray',origin='lower',aspect='auto')
plt.title('Antisymmetric mode dispersion diagram of '+str(Thickness)+' mm '+MaterialName,fontsize=20) # Rayleigh-Lamb equation amplitude in dB
plt.xlabel('Frequency (kHz)',fontsize=20)
plt.ylabel('Phase velocity (m/ms)',fontsize=20)
plt.tick_params(labelsize=15)
