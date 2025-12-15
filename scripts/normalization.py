
import numpy as np
from pathlib import Path 



#Constants:
e_0 = 8.8541878128e-12
q_e = 1.60217662e-19
k_b = 1.38064852e-23


#Box dimensions
L = 15e-3

#Dust params
Zd = 1e4
dt = 200e-6
t_max = 1
density = 1e2 #Dust particle material density
n = 6e5
lmd = 300e-6
Q = Zd * q_e
Td = 300
particle_radius = [0.001, 0.01]
m = density*np.mean(particle_radius)**3
print(m)


a = (3/(4*np.pi*n))**(1/3)
Lx = L/a
k = a/lmd


w_d = (Q*Q*n)/(4*e_0*m*a*a*a)
# t = w_d * dt /np.sqrt(3)
t = w_d * t_max


# ========================================
#   Temperature:
#
# ========================================
T = k_b * Td * 4*np.pi*e_0*a/(Q*Q)

# ========================================
#   charge:
#
#   Q = Q * (4*pi*e_0)/(k_b*T)
#
# ========================================

Qd = Q * np.sqrt((4*np.pi*e_0)/(k_b * Td))
# print("Qd=", Q)

# ========================================
#   density:
#

nd = n * a * a * a
# print("nd=", nd)
# ========================================

# coefficient for second term 
# ===========================

E_coeff = (4*np.pi*e_0*a*a*a)/ Q


# coefficient for third term 
# ===========================

g_coeff = (4*np.pi*e_0*m*a*a)/ (Q* Q)

#========================

out_str = "Normalized quantities"
out_str += "\n\n---- INPUTS -----\n"
out_str += f"\ndust radius = {particle_radius}"
out_str += f"\ndensity = {density:.2e} "
out_str += f"\nn = {n:.2e}"
out_str += f"\nlmd = {lmd:.2e}"
out_str += f"\nTemperature = {Td:.2e}"
out_str += f"\ndt = {dt:.2e}"

out_str += "\n\n ---- OUTPUTS ----\n"
out_str += f"\na = {a:.2e}"
out_str += f"\nk = {k:.2e}"
out_str += f"\ndt = {t:.2e}"
out_str += f"\nt_max = {t:.2e}"
out_str += f"\nn = {nd:.2e}"
out_str += f"\nwd = {w_d:.2e}"
out_str += f"\nQ = {Qd:.2e}"
out_str += f"\nTemperature = {T:.2e}"
out_str += f"\nLx = {Lx:.2e}"
out_str += f"\nParticle mass = {m:.2e}"
out_str += f"\nEcoeff = {E_coeff:.2e}"
out_str += f"\ngcoeff = {g_coeff:.2e}"
out_str += f"\nmass = {m:.2e}"





writefile = True
if writefile:
    out_path = Path('data') / 'normalization.txt'
    out_path.parent.mkdir(parents=True, exist_ok=True)  # ensures ./data exists
    with out_path.open('w', encoding='utf-8') as f:
        f.write(out_str)
