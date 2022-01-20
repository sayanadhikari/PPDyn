import numpy as np



#Constants:
e_0 = 8.8541878128e-12
q_e = 1.60217662e-19
k_b = 1.38064852e-23
R = 8.31446261815324 #Gas const
c = 3e8
u = 1.660539e-27

#neutral params
Tn = 300 #Temperature Kelvin
mn = 39.95 * u #atomic mass
delta = 1.44
pressure = 30 #Pa


#Dust params
Zd = 1e4
dt = 200e-6
t_max = 1
density = 1
r_d = 3.15e-6
m = density * r_d * r_d *r_d
n = 6000
lmd = 300e-6
Q = Zd * q_e
Td = 300
particle_radius = [3.25e-6,3.75e-6]

#Box dimensions
L = 70e-3

# ========================================
#   Distance:
#
#   a^2 = 1/(pi*n)
#   r = 1/a
#   k = a/lmd     screening parameter
#
# ========================================

a = np.sqrt(1/(np.pi*n))
Lx = L/a
k = a/lmd
p_r = particle_radius/a

# ========================================
#   Time:
#
#   w_d = (Q*Q*n)/(4*e_0*m*a)
#   t = w_d * t
#
# ========================================

w_d = (Q*Q*n)/(4*e_0*m*a)
t = w_d * dt
t_max_out = w_d *t_max

# ========================================
#   Temperature:
#
#   T = k_b * T * (4*pi*e_0*a)/Q^2
#
# ========================================
T = k_b * Td * 4*np.pi*e_0*a/(Q*Q)

# ========================================
#   charge:
#
#   Q = Q * (4*pi*e_0)
#
# ========================================

Qd = Q * np.sqrt((4*np.pi*e_0)/(k_b * Td))

# ========================================
#   density:
#
#   n = n * a * a
#
# ========================================

nd = n * a * a
# ========================================

# ========================================
#   Energy:
#
#   E = 4*pi*e_0*a/Q^2
#
# ========================================

K_E = 4*np.pi*e_0*a/(Q*Q)
K_m = K_E/(c*c)

# ========================================
#   Force:
#
#   F = m*a (kg m/s^2) ==> K_m * 1/(a*s*s)
#
# ========================================

K_F = K_m /(a*w_d*w_d)

#===== Neutral drag  =========


V = L*L*L
v_th_neutrals = np.sqrt((8*k_b*Tn)/(np.pi*mn))
N_neutrals = (pressure*V*6.02214076e23)/(R*Tn)
Kn_drag = delta*1.33*np.pi*N_neutrals*v_th_neutrals*mn

F_drag = Kn_drag*K_F


out_str = "Normalized quantities"
out_str += "\n\n---- INPUTS -----\n"
out_str += f"\ndust radius = {r_d}"
out_str += f"\ndensity = {density} "
out_str += f"\nn = {n}"
out_str += f"\nlmd = {lmd}"
out_str += f"\nTemperature = {Td}"
out_str += f"\ndt = {dt}"
out_str += f"\nt_max = {t_max}"
out_str += "\n\n ---- OUTPUTS ----\n"
out_str += f"\na = {a}"
out_str += f"\nk = {k}"
out_str += f"\ndt = {t}"
out_str += f"\nt_max = {t_max_out}"
out_str += f"\nn = {nd}"
out_str += f"\nwd = {w_d}"
out_str += f"\nQ = {Qd}"
out_str += f"\nTemperature = {T}"
out_str += f"\nForce coefficient = {K_F}"
out_str += f"\nMass coefficient = {K_m}"
out_str += f"\nLx = {Lx}"
out_str += f"\nParticle radius = {p_r}"
out_str += f"\nK_drag = {F_drag}"

print(out_str)