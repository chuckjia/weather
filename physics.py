from numpy import (
    exp,
    log,
)
import numpy as np


# Following are physical constants

c_p = 1004.68506  # Specific heat at constant pressure
g = 9.8
L_v = 2.53e6
L_s = 2.84e6
N_d = 200.0  # Concentration of cloud droplets
R0 = 8.31582991  # Universal gas constant
R_d = 287.058  # Gas constant for dry air
R_v = 461.5  # Gas constant for water vapor
rho_w = 1e3  # Water density

e_oo = 611.0
epsilon = 0.622  # epsilon = R_d / R_v

# Following are physical constants computed from other constants

# Relative dispersion of cloud droplet population
D_d = 0.146 - 5.964e-2 * log(N_d / 2000.0)


def rho_o_fcn(z):
    rho_oo = 1.225

    return -1e-5 * rho_oo  + rho_oo


T0 = 288.15
p0 = 101325.0

c1_p_e = g / (c_p * T0)
M_p_e = 0.02896968
c2_p_e = c_p * M_p_e / R0


def p_e_fcn(z):
    return p0 * (1 - c1_p_e * z) ** c2_p_e  # Eq.201910152104


def theta_e_fcn(z):
    theta_e_0 = 288.15

    return theta_e_0 * exp(1e-5 * z)


p_oo = 1e5
c_theta_to_T = R_d / c_p


def theta_to_T(theta, p_e):
    return theta * (p_e / p_oo) ** c_theta_to_T


def T_to_theta(T, p_e):
    return T * (p_oo / p_e) ** c_theta_to_T
