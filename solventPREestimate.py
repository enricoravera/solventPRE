import numpy as np
import scipy.constants as constants


B0 = float(input("Enter the magnetic field strength in Tesla (default 14.1): ") or 14.1)
d = float(input("Enter the distance of closest approach between the nucleus and the solvent in Angstroms (default 4.2): ") or 4.2)
d *= 1e-10 #convert to meters
S = float(input("Enter the multiplicity of the paramagnetic cosolute (default 3.5 for a gadolinium ion): ") or 3.5)
T1e = float(input("Enter the electron T1 relaxation time in seconds (default 1e-9): ") or 1e-9)
T2e = float(input("Enter the electron T2 relaxation time in seconds (default 1e-9): ") or 1e-9)
M = float(input("Enter the molar concentration of the paramagnetic cosolute in mM (default 1): ") or 1)
R2_dia = float(input("Enter the diamagnetic R2 relaxation rate in s^-1 (default 10): ") or 10)
gamma = constants.physical_constants['proton gyromag. ratio'][0] #in rad/s/T
omega_e = constants.physical_constants['electron gyromag. ratio'][0] * B0 
omega_n = gamma * B0 #in rad/s

print(f"Electron Larmor frequency: {omega_e/(2*np.pi):.2e} Hz")
print(f"Proton Larmor frequency: {omega_n/(2*np.pi):.2e} Hz")
#print(omega_e/omega_n)

k1 = (2/15) * (constants.mu_0/(4*np.pi))**2 * constants.physical_constants['Avogadro constant'][0] * (4 * np.pi /3) * gamma**2 * constants.physical_constants['Bohr magneton'][0]**2 * constants.physical_constants['electron g factor'][0]**2 * S*(S+1) * (7*T2e/(1+(omega_e*T2e)**2) + 3*T1e/(1+(omega_n*T1e)**2))
k2 = (1/15) * (constants.mu_0/(4*np.pi))**2 * constants.physical_constants['Avogadro constant'][0] * (4 * np.pi /3) * gamma**2 * constants.physical_constants['Bohr magneton'][0]**2 * constants.physical_constants['electron g factor'][0]**2 * S*(S+1) * (4*T1e + 13*T2e/(1+(omega_e*T2e)**2) + 3*T1e/(1+(omega_n*T1e)**2))

Gamma_1 = k1 * M / d**3
Gamma_2 = k2 * M / d**3

print(f"Estimated solvent PRE contribution to R1: {Gamma_1:.2f} s^-1")
print(f"Estimated solvent PRE contribution to R2: {Gamma_2:.2f} s^-1")
print(f"Optimal Tb-Ta difference for optimizing Gamma2 measurement: {1.15/(R2_dia+Gamma_2):.5f} s")
print(f"For Gamma1 sample points close to {0.5/Gamma_1:0.5f} s")
