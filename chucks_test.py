import numpy as np
import testing
import utilities as util
import time
import physics as phys

def func(x, y=0):
    a = 3
    return x + y + a

start_time = time.time()


import _model_physical as mp
mesh, _, _, _, _ = mp.param.initialize_objects()

Delta_T = 50
p_e = phys.p_e_fcn(mesh.z_ct_m)
T0 = phys.T0
init_T = T0 - (1 - p_e / phys.p0) * Delta_T
init_r_theta = phys.rho_o_fcn(mesh.z_ct_m) * phys.T_to_theta(init_T, p_e)
init_r_theta /= mesh.rho_o_m
T = phys.theta_to_T(init_r_theta, p_e)
util.print_matrix(T)



testing.show_time_used(start_time)
