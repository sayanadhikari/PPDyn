"""
Calculation of combination of constants for simulations

SI units are always used
"""
import numpy as np

#====================================================
#    Fundamental constants with convenient units
#====================================================

pi = np.pi
epsilon_0 = 8.8541878128 * 10e-12
k_b = 1.380649 * 10e-23
q_e = 1.602176634e-19 #C
u = 1.660539e-27

#====================================================
#  
#       Charge of dust particle:
#
#       Zd = a * T_e * K_zd, 
#       K_zd = 4*pi*epsilon_0*k_b/e^2
#
#====================================================


K_zd = 4*pi*epsilon_0*k_b/(q_e*q_e)# /Km such that Zd becomes unitless

if __name__ == '__main__':
    print('K_zd = ', K_zd)